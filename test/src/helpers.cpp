
#include "helpers.hpp"
#include <cstdlib>
#include <memory>
#include <random>

template <typename T> std::unique_ptr<T[]> compc_test::get_random_array(std::size_t length) {
  std::random_device r;
  std::default_random_engine e1(r());
  std::uniform_int_distribution<int> uniform_dist(1, 1000);
  std::unique_ptr<T[]> array(new T[length]);
  for (std::size_t i = 0; i < length; i++) {
    array[i] = uniform_dist(e1);
  }
  return array;
}

// template<> int16_t* get_random_array(std::size_t);
// template<> uint16_t* get_random_array<uint16_t>(std::size_t);
// template int32_t* get_random_array<int32_t>(std::size_t);
// template uint32_t* get_random_array<uint32_t>(std::size_t);
template std::unique_ptr<int64_t[]> compc_test::get_random_array<int64_t>(std::size_t length);
// template uint64_t* get_random_array<uint64_t>(std::size_t);
