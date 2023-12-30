#include "compc/elias_omega.hpp"

#include <math.h>
#include <omp.h>

#include <bitset>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <tuple>
#include <vector>

#include "compc/helpers.hpp"

template<typename T>
std::size_t compc::EliasOmega<T>::get_compressed_length(const T* array, std::size_t length)
{
  compc::ArrayPrefixSummary prefix_tuple = get_prefix_sum_array(array, length);
  return prefix_tuple.local_sums[prefix_tuple.total_chunks - 1];
}

template<typename T>
compc::ArrayPrefixSummary compc::EliasOmega<T>::get_prefix_sum_array(const T* array, std::size_t length)
{
  int local_threads = this->num_threads;
  uint32_t batch_size = this->batch_size_small;
  ;
  if (length < batch_size * static_cast<uint32_t>(local_threads))
  {  // inefficient for lenght close this
    local_threads = static_cast<int>((length + batch_size - 1) / batch_size);
  }
  else if (length >= 2 * this->batch_size_large * static_cast<uint32_t>(this->num_threads))
  {
    batch_size = this->batch_size_large;
  }
  std::size_t total_chunks = (length + batch_size - 1) / batch_size;
  std::vector<std::size_t> local_sums(total_chunks);

  bool error = false;

// manual implemenation of a omp for in order to prevent cache thrashing!
#pragma omp parallel shared(local_sums, error) num_threads(local_threads)
  {
    bool error_local = false;
    std::size_t thread_num = static_cast<std::size_t>(omp_get_thread_num());
    std::size_t num_threads_local = static_cast<std::size_t>(omp_get_num_threads());
    std::size_t start = thread_num * batch_size;
    while (true)
    {
      std::size_t l_sum = 0;
      std::size_t end = start + batch_size;
      if (start >= length)
      {
        break;
      }
      if (end > length)
      {
        end = length;
      }
      for (std::size_t i = start; i < end; i++)
      {
        T elem = array[i];
        error |= !elem;  // checking for negative inputs
        int N = hlprs::log2(static_cast<unsigned long long>(elem));
        // TODO: test for 0 and negative numbers
        while (N >= 1)
        {
          l_sum += static_cast<std::size_t>(N + 1);
          N = hlprs::log2(static_cast<unsigned long long>(N));
        }
        l_sum++;
      }
      local_sums[start / batch_size] = l_sum;
      start += num_threads_local * batch_size;
    }
#pragma omp atomic
    error |= error_local;
  }
  // final serial loop to create prefix
  // untroll --> takes almost no time, not worth it.
  std::size_t i_low = 0;
  for (std::size_t i = 1; i < total_chunks; i++)
  {
    local_sums[i] += local_sums[i_low];
    i_low = i;
  }
  return compc::ArrayPrefixSummary{ local_threads, local_sums, batch_size, total_chunks, error };  // this should use elision
}

template<typename T>
uint8_t* compc::EliasOmega<T>::compress(T* input_array, std::size_t& size)
{
  const uint64_t N = size;
  T* array;
  bool not_transformed = true;
  if (this->map_negative_numbers)
  {
    T* heap_copy_array = new T[size];
    std::memcpy(static_cast<void*>(heap_copy_array), static_cast<void*>(input_array), size * sizeof(T));
    this->transform_to_natural_numbers(heap_copy_array, size);
    array = heap_copy_array;
    not_transformed = false;
  }
  if (this->offset != 0)
  {
    if (not_transformed)
    {
      T* heap_copy_array = new T[size];
      std::memcpy(static_cast<void*>(heap_copy_array), static_cast<void*>(input_array), size * sizeof(T));
      array = heap_copy_array;
    }
    this->add_offset(array, size, this->offset);
    not_transformed = false;
  }
  if (not_transformed)
  {
    // For performance reasons we only copy the array in case map_negative_numbers is true;
    array = input_array;
  }

  ArrayPrefixSummary prefix_tuple = this->get_prefix_sum_array(array, N);  // in bits
  if (prefix_tuple.error)
  {
    if (!not_transformed)
    {
      delete[] array;
    }
    return nullptr;
  }

  int local_threads = prefix_tuple.local_threads;
  std::vector<std::size_t> prefix_array = prefix_tuple.local_sums;
  uint32_t batch_size = prefix_tuple.batch_size;
  std::size_t total_chunks = prefix_tuple.total_chunks;
  const uint64_t compressed_length = prefix_array[total_chunks - 1];
  const uint64_t compressed_bytes = (compressed_length + 7) / 8;  // getting the number of bytes (ceil)
  // zero initialize, otherwise there are problems at the edges of the batches
  uint8_t* compressed = new uint8_t[compressed_bytes]();

#pragma omp parallel shared(compressed, prefix_array) firstprivate(N) num_threads(local_threads)
  {
    std::size_t start_bit;
    std::size_t start_index;
#pragma omp for schedule(dynamic, batch_size)
    for (uint32_t round = 0; round < total_chunks; round++)
    {
      uint8_t current_byte = 0;
      if (round == 0)
      {
        start_bit = 0;
        start_index = 0;
      }
      else
      {
        start_bit = prefix_array[round - 1];
        start_index = round * batch_size;
      }
      std::size_t end_bit = prefix_array[round];
      std::size_t end_index = start_index + batch_size;
      if (end_index > N)
      {
        end_index = N;
      }
      std::size_t start_byte = start_bit / 8;
      std::size_t end_byte = end_bit / 8;
      std::size_t index = start_byte;  // index for the byte array
      int bits_left = 8 - (static_cast<int>(start_bit) - static_cast<int>(start_byte) * 8);
      for (std::size_t i = start_index; i < end_index; i++)
      {
        T local_N = array[i];
        // unrolling recursive definition of omega coding
        std::vector<T> v;
        v.push_back(0);

        while (local_N > 1)
        {
          v.push_back(local_N);
          T N_binary_length = static_cast<T>(hlprs::log2(static_cast<unsigned long long>(local_N))) + 1;
          local_N = N_binary_length - 1;
        }

        T local_binary_length = 0;
        T old_local_value = 1;
        for (auto it = v.rbegin(); it != v.rend(); ++it)
        {
          T local_value = *it;
          if (local_value == 0)
          {
            local_binary_length = 1;
          }
          else
          {
            local_binary_length = old_local_value + 1;
          }
          old_local_value = local_value;

          while (local_binary_length)
          {
            uint8_t mask = 255u;
            mask = mask >> (8 - bits_left);
            if (bits_left > 0 && local_binary_length >= static_cast<T>(bits_left))
            {
              local_binary_length = local_binary_length - static_cast<T>(bits_left);
              current_byte = current_byte | static_cast<uint8_t>((local_value >> local_binary_length) & mask);
              bits_left = 0;
              if (index == start_byte || index == end_byte)
              {
#pragma omp atomic
                compressed[index] |= current_byte;
              }
              else
              {
                compressed[index] = current_byte;
              }
              index++;
              current_byte = 0;
              bits_left = 8;
            }
            else if (bits_left > 0 && local_binary_length < static_cast<T>(bits_left))
            {
              current_byte =
                  current_byte |
                  static_cast<uint8_t>((local_value << (bits_left - static_cast<uint8_t>(local_binary_length))) & mask);
              bits_left -= static_cast<uint8_t>(local_binary_length);
              local_binary_length = 0;
            }
          }
        }
      }
      if (bits_left < 8)
      {
#pragma omp atomic
        compressed[index] |= current_byte;
      }
    }
  }
  std::cout << "fin" << std::endl;
  if (this->map_negative_numbers)
  {
    delete[] array;
  }
  size = compressed_bytes;
  return compressed;
}

template<typename T>
T* compc::EliasOmega<T>::decompress(const uint8_t* array, std::size_t binary_length, std::size_t array_length)
{
  T* uncomp = new T[array_length];
  std::size_t index = 0;
  T current_decoded_number = 0;
  std::size_t binary_index = 0;
  uint8_t current_byte = array[binary_index];
  uint8_t bits_left = 8;
  binary_length -= 1;
  bool read_new_number = true;
  std::size_t to_read_left = 0;
  while (binary_index < binary_length || index < array_length)
  {
    if (bits_left == 0)
    {
      binary_index++;
      current_byte = array[binary_index];
      bits_left = 8;
    }
    uint8_t cur_copy = current_byte;
    current_byte = current_byte << (8 - bits_left);
    if (read_new_number)
    {
      bool state = (current_byte >> 7);
      if (!state)
      {
        uncomp[index] = 1;
        index++;
        bits_left--;
        current_byte = cur_copy;
        continue;
      }
      read_new_number = false;
      to_read_left = 2;
    }
    if (to_read_left == 0)
    {
      bool state = (current_byte >> 7);
      if (!state)
      {
        uncomp[index] = current_decoded_number;
        index++;
        bits_left--;
        current_byte = cur_copy;
        read_new_number = true;
        current_decoded_number = 0;
        continue;
      }
      else
      {
        to_read_left = static_cast<std::size_t>(current_decoded_number + 1);
        current_decoded_number = 0;
      }
    }
    current_byte = cur_copy;
    while (to_read_left && bits_left)
    {
      uint8_t mask = 255u >> (8 - bits_left);
      T curT = static_cast<T>(current_byte & mask);
      uint8_t bits_to_process =
          (to_read_left < static_cast<std::size_t>(bits_left)) ? static_cast<uint8_t>(to_read_left) : bits_left;
      bits_left -= bits_to_process;
      to_read_left -= bits_to_process;
      current_decoded_number = current_decoded_number | static_cast<T>((curT << to_read_left) >> bits_left);
    }
  }
  if (this->offset != 0)
  {
    this->add_offset(uncomp, array_length, -this->offset);
  }
  if (this->map_negative_numbers)
  {
    this->transform_to_natural_numbers_reverse(uncomp, array_length);
  }
  return uncomp;
}

template class compc::EliasOmega<int16_t>;
template class compc::EliasOmega<uint16_t>;
template class compc::EliasOmega<int32_t>;
template class compc::EliasOmega<uint32_t>;
template class compc::EliasOmega<int64_t>;
template class compc::EliasOmega<uint64_t>;