#ifndef COMPC_ELIAS_DELTA_H_
#define COMPC_ELIAS_DELTA_H_
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>

#include "compintc/elias_base.hpp"
namespace compc {

template <typename T> class EliasDelta : public EliasBase<T> {
public:
  uint32_t batch_size_small{50};
  uint32_t batch_size_large{1000};
  EliasDelta() = default;
  EliasDelta(T zero_offset, bool map_negative_numbers_to_positive)
      : EliasBase<T>(zero_offset, map_negative_numbers_to_positive){};
  EliasDelta(T zero_offset, bool map_negative_numbers_to_positive, uint32_t batch_size_small_p,
             uint32_t batch_size_large_p)
      : EliasBase<T>(zero_offset, map_negative_numbers_to_positive), batch_size_small(batch_size_small_p),
        batch_size_large(batch_size_large_p){};
  ~EliasDelta() = default;
  std::unique_ptr<uint8_t[]> compress(const T*, std::size_t&) override;
  std::unique_ptr<T[]> decompress(const uint8_t*, std::size_t, std::size_t) override;
  std::size_t get_compressed_length(const T*, std::size_t) override;
  ArrayPrefixSummary get_prefix_sum_array(const T*, std::size_t) override;
  // copy constructor
  EliasDelta(EliasDelta& other)
      : EliasDelta(other.offset, other.map_negative_numbers, other.batch_size_small, other.batch_size_large){};
  // move constructor
  EliasDelta(EliasDelta&& other) noexcept
      : EliasBase<T>(other), batch_size_small(std::exchange(batch_size_small, 0)),
        batch_size_large(std::exchange(batch_size_large, 0)){};
  EliasDelta& operator=(EliasDelta other) {
    this->num_threads = other.num_threads;
    this->offset = other.offset;
    this->map_negative_numbers = other.map_negative_numbers;
    batch_size_small = other.batch_size_small;
    batch_size_large = other.batch_size_large;
    return *this;
  };
  EliasDelta& operator=(EliasDelta&& other) noexcept {
    this->num_threads = std::move(other.num_threads);
    this->offset = std::move(other.offset);
    this->map_negative_numbers = std::move(other.map_negative_numbers);
    batch_size_small = std::move(other.batch_size_small);
    batch_size_large = std::move(other.batch_size_large);
    return *this;
  };
};
} // namespace compc

#endif // COMPC_ELIAS_DELTA_H_
