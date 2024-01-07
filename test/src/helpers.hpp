#ifndef COMPC_TESTS_HELPERS_H_
#define COMPC_TESTS_HELPERS_H_
#include <cstdint>
#include <cstdlib>
#include <memory>
namespace compc_test {
template <typename T> std::unique_ptr<T[]> get_random_array(std::size_t length);
}

#endif // COMPC_TESTS_HELPERS_H_
