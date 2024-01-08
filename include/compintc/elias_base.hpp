#ifndef COMPC_ELIAS_BASE_H_
#define COMPC_ELIAS_BASE_H_
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "compintc/compressor.hpp"
#include "compintc/helpers.hpp"
namespace compc {

struct ArrayPrefixSummary {
  int local_threads = 0;
  uint32_t batch_size = 0;
  std::vector<std::size_t> local_sums{};
  std::size_t total_chunks = 0;
  bool error = false;
};

template <typename T> class EliasBase : public Compressor<T> {
public:
  T offset{0};
  bool map_negative_numbers{false};
  EliasBase() = default;
  explicit EliasBase(T zero_offset) : offset(zero_offset){};
  EliasBase(T zero_offset, bool map_negative_numbers_to_positive)
      : offset(zero_offset), map_negative_numbers(map_negative_numbers_to_positive){};
  virtual ~EliasBase() = default;
  virtual ArrayPrefixSummary get_prefix_sum_array(const T*, std::size_t) = 0;
  // copy constructor
  EliasBase(EliasBase& other)
      : Compressor<T>(other), offset(other.offset), map_negative_numbers(other.map_negative_numbers){};
  // move constructor
  EliasBase(EliasBase&& other) noexcept // move constructor
      : Compressor<T>(other), offset(std::exchange(other.offset, 0)),
        map_negative_numbers(std::exchange(other.map_negative_numbers, false)){};
  // copy operator
  EliasBase& operator=(const EliasBase& other) = default;
  EliasBase& operator=(EliasBase&& other) noexcept = default;

protected:
  std::unique_ptr<T[]> transform_array_inputs(const T* input_array, std::size_t& size) {
    std::unique_ptr<T[]> heap_copy_array = nullptr; // TODO change to make_unique_for_overwrite
    if (this->map_negative_numbers || this->offset != 0) {
      heap_copy_array = std::move(std::unique_ptr<T[]>(new T[size]));
    }
    bool not_transformed = true;
    if (this->map_negative_numbers) {
      std::memcpy(static_cast<void*>(heap_copy_array.get()), static_cast<const void*>(input_array), size * sizeof(T));
      this->transform_to_natural_numbers(heap_copy_array.get(), size);
      if (this->offset != 0) {
        this->add_offset(heap_copy_array.get(), size, this->offset);
      }
      not_transformed = false;
    }
    if (not_transformed && this->offset != 0) {
      std::memcpy(static_cast<void*>(heap_copy_array.get()), static_cast<const void*>(input_array), size * sizeof(T));
      this->add_offset(heap_copy_array.get(), size, this->offset);
    }
    return heap_copy_array;
  }
};
} // namespace compc

#endif // COMPC_ELIAS_BASE_H_
