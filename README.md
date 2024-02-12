# CompIntC

A high-performance C++ library for compressing integers with variable-length Elias compression algorithms, which are examples of [universal codes](https://en.wikipedia.org/wiki/Universal_code_(data_compression)). The library was designed for compressing indices of weight matrices in decentralised machine learning work loads. Hence, it supports multithreading for faster encoding of indices to speed up the distribution of the data to neighbouring nodes.


## Supported Algorithms


The library implements the Elias [gamma](https://en.wikipedia.org/wiki/Elias_gamma_coding), [delta](https://en.wikipedia.org/wiki/Elias_delta_coding), and [omega](https://en.wikipedia.org/wiki/Elias_omega_coding) algorithms.


All three are in the `compc` namespace and can be used like this (Replace `gamma` with your respective choices):

```
#include "compc/elias_gamma.hpp"
compc::EliasGamma<long> elias;
std::size_t size = 10;
long input[10] =  {1, 3, 2000, 2, 50, 1,25345, 11, 10000000, 1};
uint8_t* comp = elias.compress(input, size);
long* output = elias.decompress(comp, size, 10);
```

## Documentation


All three algorithms are based on a common base class `EliasBase\<typename T\>` and support off-sets and mappings for negative numbers, which are needed to compress all numbers instead of just positive integers.


The base class has the following constructors constructors:

```
   EliasBase()
   EliasBase(T zero_offset)
   EliasBase(T zero_offset, bool
```

with:
- zero_offset: Offset from zero. Is applied before the input array is compressed.
- map_negative_numbers_to_positive: If set negative numbers are mapped to positive numbers before compression.

Its subclasses can be initiated for the following types:
- int16_t
- uint16_t
- int32_t
- uint32_t
- int64_t
- uint64_t


Furthermore, there are two methods for compression and decompression:

```
uint8_t* compress(const T* input_array, std::size_t& size)
```
- input_array: The array to be compressed.
- size: The size of the input_array. It gets overwritten by the size of the produced output.

```
T* decompress(const uint8_t* array, std::size_t binary_length, std::size_t array_length)
```
- array: Compressed array
- binary_length: Size of the compressed array
- array_length: Size of the output array.


## Multi-threading
The compress function is parallelized with OpenMP. You can set the number of threads by setting the `OMP_NUM_THREADS` environment variable, e.g.,
```
export OMP_NUM_THREADS=10
```

Alternatively, you can overwrite this global choice by setting:
```
compc::EliasGamma<long> elias;
elias->num_threads = 5;
```

## Bindings

There exist Python bindings for the library. See our sister project [ComIntPy](https://github.com/JeffWigger/compintpy).
