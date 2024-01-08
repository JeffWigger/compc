#include "compintc/elias_gamma.hpp"

#include <cmath>
#include <omp.h>

#include <bitset>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

#include "compintc/helpers.hpp"

template <typename T> std::size_t compc::EliasGamma<T>::get_compressed_length(const T* array, std::size_t length) {
  compc::ArrayPrefixSummary prefix_tuple = get_prefix_sum_array(array, length);
  return prefix_tuple.local_sums[prefix_tuple.total_chunks - 1];
}

template <typename T>
compc::ArrayPrefixSummary compc::EliasGamma<T>::get_prefix_sum_array(const T* array, std::size_t length) {
  int local_threads = this->num_threads;
  uint32_t batch_size = this->batch_size_small;
  // inefficient for lenght close this
  if (length < static_cast<std::size_t>(batch_size) * static_cast<std::size_t>(local_threads)) {
    local_threads = static_cast<int>((length + batch_size - 1) / batch_size);
  } else if (length >= 2 * this->batch_size_large * static_cast<uint32_t>(this->num_threads)) {
    batch_size = this->batch_size_large;
  }
  std::size_t total_chunks = (length + batch_size - 1) / batch_size;
  std::vector<std::size_t> local_sums(total_chunks);

  bool error = false;

// manual implemenation of a omp for in order to prevent cache thrashing!
#pragma omp parallel default(none) shared(local_sums, error, array) firstprivate(batch_size, length)                   \
    num_threads(local_threads)
  {
    bool error_local = false;
    auto thread_num = static_cast<std::size_t>(omp_get_thread_num());
    auto num_threads_local = static_cast<std::size_t>(omp_get_num_threads());
    std::size_t start = thread_num * batch_size;
    while (true) {
      std::size_t l_sum = 0;
      std::size_t end = start + batch_size;
      if (start >= length) {
        break;
      }
      if (end > length) {
        end = length;
      }
      //#pragma omp unroll partial(4)
      for (std::size_t i = start; i < end; i++) {
        T elem = array[i];
        error |= !elem; // checking for negative inputs
        // 2*N + 1
        l_sum += (static_cast<uint64_t>(hlprs::log2(static_cast<unsigned long long>(elem))) << 1U) + 1;
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
  //#pragma omp unroll partial(4)
  for (std::size_t i = 1; i < total_chunks; i++) {
    local_sums[i] += local_sums[i_low];
    i_low = i;
  }
  return compc::ArrayPrefixSummary{local_threads, batch_size, local_sums, total_chunks,
                                   error}; // this should use elision
}

template <typename T>
std::unique_ptr<uint8_t[]> compc::EliasGamma<T>::compress(const T* input_array, std::size_t& size) {
  const uint64_t N = size;
  const T* array = nullptr;
  std::unique_ptr<T[]> heap_copy_array; // TODO change to make_unique_for_overwrite
  if (this->map_negative_numbers || this->offset != 0) {
    heap_copy_array = this->transform_array_inputs(input_array, size);
    array = heap_copy_array.get();
  } else {
    array = input_array;
  }
  ArrayPrefixSummary prefix_tuple = this->get_prefix_sum_array(array, N); // in bits
  if (prefix_tuple.error) {
    return nullptr;
  }

  int local_threads = prefix_tuple.local_threads;
  std::vector<std::size_t> prefix_array = prefix_tuple.local_sums;
  uint32_t batch_size = prefix_tuple.batch_size;
  std::size_t total_chunks = prefix_tuple.total_chunks;
  const uint64_t compressed_length = prefix_array[total_chunks - 1];
  const uint64_t compressed_bytes = (compressed_length + 7) / 8; // getting the number of bytes (ceil)
  // zero initialize, otherwise there are problems at the edges of the batches
  std::unique_ptr<uint8_t[]> compressed = std::make_unique<uint8_t[]>(compressed_bytes);

#pragma omp parallel default(none) shared(compressed, prefix_array, array) firstprivate(N, total_chunks, batch_size)   \
    num_threads(local_threads)
  {
    std::size_t start_bit = 0;
    std::size_t start_index = 0;
#pragma omp for schedule(dynamic, batch_size)
    for (uint32_t round = 0; round < total_chunks; round++) {
      uint8_t current_byte = 0;
      if (round == 0) {
        start_bit = 0;
        start_index = 0;
      } else {
        start_bit = prefix_array[round - 1];
        start_index = static_cast<std::size_t>(round) * static_cast<std::size_t>(batch_size);
      }
      std::size_t end_bit = prefix_array[round];
      std::size_t end_index = start_index + batch_size;
      if (end_index > N) {
        end_index = N;
      }
      std::size_t start_byte = start_bit / 8;
      std::size_t end_byte = end_bit / 8;
      std::size_t index = start_byte; // index for the byte array
      uint bits_left = 8 - (static_cast<uint>(start_bit) - static_cast<uint>(start_byte) * 8);
      for (std::size_t i = start_index; i < end_index; i++) {
        T value = array[i];
        uint length_prefix_part = static_cast<uint>(hlprs::log2(static_cast<unsigned long long>(value)));
        uint length_binary_part = length_prefix_part + 1;

        // Part 1: writing the prefix 0s
        while (length_prefix_part > 0) {
          if (bits_left > length_prefix_part) {
            bits_left -= length_prefix_part;
            length_prefix_part = 0;
          } else {
            if (index == start_byte) {
#pragma omp atomic
              compressed[index] = compressed[index] | current_byte;
            } else {
              compressed[index] = compressed[index] | current_byte;
            }
            index++;
            length_prefix_part = length_prefix_part - bits_left;
            current_byte = 0;
            bits_left = 8;
          }
        }
        // Part 2: writing the number in binary
        while (length_binary_part > 0) {
          uint8_t mask = 255U;
          mask = mask >> (8U - bits_left);
          if (bits_left > 0 && length_binary_part >= bits_left) {
            length_binary_part = length_binary_part - bits_left;
            current_byte = current_byte | static_cast<uint8_t>((value >> length_binary_part) & mask);
            bits_left = 0;
            if (index == start_byte || index == end_byte) {
#pragma omp atomic
              compressed[index] |= current_byte;
            } else {
              compressed[index] = current_byte;
            }
            index++;
            current_byte = 0;
            bits_left = 8;
          } else if (bits_left > 0 && length_binary_part < bits_left) {
            current_byte = current_byte | static_cast<uint8_t>((value << (bits_left - length_binary_part)) & mask);
            bits_left -= length_binary_part;
            length_binary_part = 0;
          }
        }
      }
      if (bits_left < 8) {
#pragma omp atomic
        compressed[index] |= current_byte;
      }
    }
  }
  size = compressed_bytes;
  return compressed;
}

template <typename T>
std::unique_ptr<T[]> compc::EliasGamma<T>::decompress(const uint8_t* array, std::size_t binary_length,
                                                      std::size_t array_length) {
  std::unique_ptr<T[]> uncomp(new T[array_length]);
  std::size_t index = 0;
  T current_decoded_number = 0;
  uint length_binary_part = 0;
  bool reading_prefix_zeros = true;
  std::size_t binary_index = 0;
  uint8_t current_byte = array[binary_index];
  uint8_t bits_left = 8;
  binary_length -= 1;
  while (binary_index < binary_length || index < array_length) {
    if (bits_left == 0) {
      binary_index++;
      current_byte = array[binary_index];
      bits_left = 8;
    }
    uint8_t cur_copy = current_byte;
    current_byte = current_byte << (8U - bits_left);
    // TODO: process more than one bit at a time.
    while (reading_prefix_zeros && bits_left > 0) {
      bool state = !(current_byte >> 7U); // is either 0, or 1
      current_byte = current_byte << state;
      bits_left -= state;
      reading_prefix_zeros = state;
      length_binary_part++;
    }
    current_byte = cur_copy;

    while (!reading_prefix_zeros && bits_left > 0) {
      uint8_t mask = 255U >> (8U - bits_left);
      T curT = static_cast<T>(current_byte & mask);
      bool state = (length_binary_part >= bits_left);
      uint8_t bits_to_process = (state) ? bits_left : static_cast<uint8_t>(length_binary_part);
      bits_left -= bits_to_process;
      length_binary_part -= bits_to_process;
      current_decoded_number = current_decoded_number | static_cast<T>((curT << length_binary_part) >> bits_left);
      reading_prefix_zeros = !length_binary_part;
      if (reading_prefix_zeros) {
        uncomp[index] = current_decoded_number;
        index++;
        current_decoded_number = 0U;
      }
    }
  }
  if (this->offset != 0) {
    this->add_offset(uncomp.get(), array_length, -this->offset);
  }
  if (this->map_negative_numbers) {
    this->transform_to_natural_numbers_reverse(uncomp.get(), array_length);
  }
  return uncomp;
}

template class compc::EliasGamma<int16_t>;
template class compc::EliasGamma<uint16_t>;
template class compc::EliasGamma<int32_t>;
template class compc::EliasGamma<uint32_t>;
template class compc::EliasGamma<int64_t>;
template class compc::EliasGamma<uint64_t>;
