#ifndef COMPC_HELPERS_H_
#define COMPC_HELPERS_H_
#include <cstdint>

namespace hlprs {
inline int log2(unsigned long long x) {
  // 0 is not a valid input
  // https://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c
  // https://stackoverflow.com/questions/70015481/how-do-i-count-leading-zeros-on-both-mac-m1-and-x86-64
  return (63 - __builtin_clzll(x)); // todo: replace 63 with sizeof(long long)*8 - 1
  // in c++20 we could use std::bit_width(index) - 1
}

} // namespace hlprs
#endif // COMPC_HELPERS_H_
