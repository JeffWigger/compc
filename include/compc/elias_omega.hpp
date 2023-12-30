#ifndef COMPC_ELIAS_OMEGA_H_
#define COMPC_ELIAS_OMEGA_H_
#include <cstdint>
#include <tuple>

#include "compc/elias_base.hpp"
namespace compc
{

  template<typename T>
  class EliasOmega : public EliasBase<T>
  {
   public:
    EliasOmega() = default;
    EliasOmega(T zero_offset, bool map_negative_numbers_to_positive)
        : EliasBase<T>(zero_offset, map_negative_numbers_to_positive){};
    ~EliasOmega() = default;
    uint8_t* compress(T*, std::size_t&) override;
    T* decompress(const uint8_t*, std::size_t, std::size_t) override;
    std::size_t get_compressed_length(const T*, std::size_t) override;
    ArrayPrefixSummary get_prefix_sum_array(const T*, std::size_t) override;

   private:
    uint32_t batch_size_small{ 50 };
    uint32_t batch_size_large{ 1000 };
  };
}  // namespace compc

#endif  // COMPC_ELIAS_DELTA_H_
