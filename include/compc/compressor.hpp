#ifndef COMPC_COMPRESSOR_H_
#define COMPC_COMPRESSOR_H_
#include <cstdlib>
#include <omp.h>

#include <cstdint>
#include <memory>
namespace compc {

template <typename T> class Compressor {
public:
  int num_threads{1};
  Compressor() {
    char* num_threads_char = std::getenv("OMP_NUM_THREADS");
    if (num_threads_char != nullptr) {
      int omp_threads = static_cast<int>(std::strtol(num_threads_char, nullptr, 10));
      this->num_threads = std::max(omp_threads, 1);
    } else {
      this->num_threads = 1;
    }
    omp_set_num_threads(this->num_threads);
  };
  explicit Compressor(int number_of_threads) : num_threads(number_of_threads) {
    omp_set_num_threads(this->num_threads);
  };
  virtual ~Compressor() = default;
  virtual std::unique_ptr<uint8_t[]> compress(const T*, std::size_t&) = 0;
  virtual std::unique_ptr<T[]> decompress(const uint8_t*, std::size_t, std::size_t) = 0;
  virtual std::size_t get_compressed_length(const T*, std::size_t) = 0;
  // copy cunstructor
  Compressor(Compressor& other) : Compressor(other.num_threads){};
  // move cunstructor
  Compressor(Compressor&& other) noexcept // move constructor
      : num_threads(std::exchange(other.num_threads, 0)){};
  // copy operator
  Compressor& operator=(const Compressor& other) = default;
  Compressor& operator=(Compressor&& other) noexcept = default;

  void transform_to_natural_numbers(T* array, const std::size_t& size)
  /*
    array: array of numbers to transform
    size: size of the array

    The array input is transformed back to numbers in situ.
    Positive number a --> 2 * a
    Negative number b --> -2 * b - 1

  */
  {
    int local_threads = this->num_threads;
    if (size < static_cast<std::size_t>(local_threads)) {
      local_threads = 1;
    }
#pragma omp parallel for schedule(guided, 1000) default(none) shared(array) firstprivate(size)                         \
    num_threads(local_threads)
    // #pragma omp unroll partial(4)
    for (std::size_t i = 0; i < size; i++) {
      T at_i = array[i];
      T bi = (at_i < 0);
      array[i] = static_cast<T>(at_i * (2 - 4 * bi)) - bi;
    }
  }

  void transform_to_natural_numbers_reverse(T* array, const std::size_t& size)
  /*
    array: array of numbers to transform
    size: size of the array

    The array input is transformed back to numbers in situ.
    Even number a --> a / 2
    Odd number b --> -(b + 1) / 2

  */
  {
    int local_threads = this->num_threads;
    if (size < static_cast<std::size_t>(local_threads)) {
      local_threads = 1;
    }
#pragma omp parallel for schedule(guided, 1000) default(none) shared(array) firstprivate(size)                         \
    num_threads(local_threads)
    // #pragma omp unroll partial(4)
    for (std::size_t i = 0; i < size; i++) {
      T at_i = array[i];
      T bi = (at_i % 2);
      array[i] = ++at_i / static_cast<T>((2 - 4 * bi));
    }
  }

  void add_offset(T* array, const std::size_t& size, T offset) {
    int local_threads = this->num_threads;
    if (size < static_cast<std::size_t>(local_threads)) {
      local_threads = 1;
    }
#pragma omp parallel for schedule(guided, 1000) default(none) shared(array) firstprivate(size, offset)                 \
    num_threads(local_threads)
    // #pragma omp unroll partial(4)
    for (std::size_t i = 0; i < size; i++) {
      array[i] = array[i] + offset;
    }
  }
};
} // namespace compc

#endif // COMPC_COMPRESSOR_H_
