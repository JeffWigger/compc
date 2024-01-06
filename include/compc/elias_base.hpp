#ifndef COMPC_ELIAS_BASE_H_
#define COMPC_ELIAS_BASE_H_
#include <cmath>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "compc/compressor.hpp"
#include "compc/helpers.hpp"
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
  EliasBase& operator=(EliasBase other) {
    *this = EliasBase(other);
    return *this;
  };
  EliasBase& operator=(EliasBase&& other) noexcept {
    Compressor<T>::operator=(std::forward(other));
    offset = std::move(other.offset);
    map_negative_numbers = std::move(other.map_negative_numbers);
    return *this;
  };
};
} // namespace compc

#endif // COMPC_ELIAS_BASE_H_
