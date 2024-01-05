#include "compc/elias_omega.hpp"
#include "helpers.hpp"
#include <bitset>
#include <chrono>
#include <cstdint>
#include <gtest/gtest.h>
using namespace std::chrono;

TEST(Elias_Omega_Unit_Delta_GetCompressedLength, CheckValues) {
  std::size_t size = 7;
  long input[7] = {1, 2, 3, 4, 5, 10, 17};
  compc::EliasOmega<long> elias;
  std::size_t length = elias.get_compressed_length(input, size);
  std::cout << "compressed length: " << length << std::endl;
  ASSERT_EQ(length, 37); // comparing values
}

TEST(Elias_Omega_Unit_Delta_GetCompressedLength_10, CheckValues) {
  std::size_t size = 11;
  long input[11] = {1, 2, 3, 4, 5, 10, 17, 100, 10000, 100000, 1000000};
  compc::EliasOmega<long> elias;
  std::size_t length = elias.get_compressed_length(input, size);
  std::cout << "compressed length: " << length << std::endl;
  ASSERT_EQ(length, 130); // comparing values
}

TEST(Elias_Omega_Unit_Encode, CheckValues) {
  std::size_t size = 5;
  long input[5] = {1, 2, 5, 10, 17};
  uint8_t output[4] = {74, 186, 82, 32};
  compc::EliasOmega<long> elias;
  uint8_t* comp = elias.compress(input, size);
  for (std::size_t i = 0; i < size; i++) {
    std::cout << std::bitset<8>(comp[i]) << std::endl;
    ASSERT_EQ(output[i], comp[i]);
  }
  delete[] comp;
}

TEST(Elias_Omega_Unit_Decode, CheckValues) {
  std::size_t size = 5;
  uint8_t input[4] = {74, 186, 82, 32};
  long output[5] = {1, 2, 5, 10, 17};
  compc::EliasOmega<long> elias;
  long* res = elias.decompress(input, 4, size);
  for (std::size_t i = 0; i < size; i++) {
    std::cout << res[i] << std::endl;
    ASSERT_EQ(output[i], res[i]);
  }

  delete[] res;
}

TEST(Elias_Omega_DecompCompEQTestLong, CheckValues) {
  std::size_t size = 10;
  long input[10] = {1, 3, 2000, 2, 50, 1, 25345, 11, 10000000, 1};
  compc::EliasOmega<long> elias;
  std::cout << "size: " << sizeof(input) << std::endl;
  uint8_t* comp = elias.compress(input, size);

  std::cout << "binary size: " << size << std::endl;
  for (std::size_t i = 0; i < size; i++)
    std::cout << std::bitset<8>(comp[i]) << std::endl;
  long* output = elias.decompress(comp, size, 10);
  // ASSERT_EQ(size, 21);
  for (int i = 0; i < 10; i++) {
    ASSERT_EQ(output[i], input[i]); // comparing values
  }
  delete[] comp;
  delete[] output;
}

TEST(Elias_Omega_SpeedTestLong, CheckValues) {
  std::size_t len = 500000;
  std::size_t len_copy = len;
  long* random_array = compc_test::get_random_array<long>(len);

  auto start = high_resolution_clock::now();
  compc::EliasOmega<long> elias;
  uint8_t* comp = elias.compress(random_array, len);
  std::cout << "binary size: " << len << std::endl;
  auto mid = high_resolution_clock::now();
  long* output = elias.decompress(comp, len, len_copy);
  // std::cout << "decompress" << std::endl;

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  auto mid_duration = duration_cast<microseconds>(mid - start);
  auto decomp_duration = duration_cast<microseconds>(stop - mid);
  std::cout << "Total: " << duration.count() << " compression: " << mid_duration.count()
            << " decompression: " << decomp_duration.count() << std::endl;
  for (std::size_t i = 0; i < len_copy; i++) {
    // std::cout << i << std::endl;
    ASSERT_EQ(output[i], random_array[i]); // comparing values
  }
  delete[] comp;
  delete[] output;
  delete[] random_array;
}

TEST(Elias_Omega_DecompCompEQTestInt, CheckValues) {
  std::size_t size = 10;
  int input[10] = {1, 3, 2000, 2, 50, 1, 25345, 11, 1000000, 1};
  compc::EliasOmega<int> elias;
  std::cout << "size: " << sizeof(input) << std::endl;
  uint8_t* comp = elias.compress(input, size);

  int* output = elias.decompress(comp, size, 10);
  // ASSERT_EQ(size, 21);
  for (int i = 0; i < 10; i++) {
    ASSERT_EQ(output[i], input[i]); // comparing values
  }
  delete[] comp;
  delete[] output;
}

TEST(Elias_Omega_DecompCompEQTestShort, CheckValues) {
  std::size_t size = 10;
  short input[10] = {1, 3, 2000, 2, 50, 1, 25345, 11, 10000, 1};
  compc::EliasOmega<short> elias;
  std::cout << "size: " << sizeof(input) << std::endl;
  uint8_t* comp = elias.compress(input, size);

  short* output = elias.decompress(comp, size, 10);
  // ASSERT_EQ(size, 21);
  for (int i = 0; i < 10; i++) {
    ASSERT_EQ(output[i], input[i]); // comparing values
  }
  delete[] comp;
  delete[] output;
}

TEST(Helpers_transform_to_natural_numbers, CheckValues) {
  std::size_t size = 8;
  long input[8] = {1, -1, -5, 5, -100, 100, 10000, -10000};
  long input_gt[8] = {1, -1, -5, 5, -100, 100, 10000, -10000};
  long results_gt[8] = {2, 1, 9, 10, 199, 200, 20000, 19999};
  compc::EliasOmega<long> elias;
  elias.transform_to_natural_numbers(input, size);
  for (int i = 0; i < 8; i++) {
    ASSERT_EQ(input[i], results_gt[i]); // comparing values
  }
  elias.transform_to_natural_numbers_reverse(input, size);
  for (int i = 0; i < 8; i++) {
    ASSERT_EQ(input[i], input_gt[i]); // comparing values
  }
}

TEST(Helpers_transform_to_natural_numbers_switch_11, CheckValues) {
  std::size_t size = 11;
  long input[11] = {1, -1, -5, 5, -100, 100, 10000, -10000, -5, 5, -100};
  long input_gt[11] = {1, -1, -5, 5, -100, 100, 10000, -10000, -5, 5, -100};
  long results_gt[11] = {2, 1, 9, 10, 199, 200, 20000, 19999, 9, 10, 199};
  compc::EliasOmega<long> elias;
  elias.transform_to_natural_numbers(input, size);
  for (int i = 0; i < 11; i++) {
    ASSERT_EQ(input[i], results_gt[i]); // comparing values
  }
  elias.transform_to_natural_numbers_reverse(input, size);
  for (int i = 0; i < 11; i++) {
    ASSERT_EQ(input[i], input_gt[i]); // comparing values
  }
}

TEST(Helpers_transform_to_natural_numbers_switch_9, CheckValues) {
#define SIZE 9
  std::size_t size = SIZE;
  long input[SIZE] = {1, -1, -5, 5, -100, 100, 10000, -10000, -6};
  long input_gt[SIZE] = {1, -1, -5, 5, -100, 100, 10000, -10000, -6};
  long results_gt[SIZE] = {2, 1, 9, 10, 199, 200, 20000, 19999, 11};
  compc::EliasOmega<long> elias;
  elias.transform_to_natural_numbers(input, size);
  for (int i = 0; i < SIZE; i++) {
    ASSERT_EQ(input[i], results_gt[i]); // comparing values
  }
  elias.transform_to_natural_numbers_reverse(input, size);
  for (int i = 0; i < SIZE; i++) {
    ASSERT_EQ(input[i], input_gt[i]); // comparing values
  }
}

TEST(Helpers_transform_to_natural_numbers_SpeedTestLong, CheckValues) {
  std::size_t len = 500000;
  long* random_array = compc_test::get_random_array<long>(len);

  auto start = high_resolution_clock::now();
  compc::EliasOmega<long> elias;
  elias.transform_to_natural_numbers(random_array, len);
  auto mid = high_resolution_clock::now();
  elias.transform_to_natural_numbers_reverse(random_array, len);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  auto mid_duration = duration_cast<microseconds>(mid - start);
  auto decomp_duration = duration_cast<microseconds>(stop - mid);
  std::cout << "Total: " << duration.count() << " forward: " << mid_duration.count()
            << " backward: " << decomp_duration.count() << std::endl;
  delete[] random_array;
}

TEST(Helpers_add_offset_SpeedTestLong, CheckValues) {
  std::size_t len = 500000;
  long* random_array = compc_test::get_random_array<long>(len);

  auto start = high_resolution_clock::now();
  compc::EliasOmega<long> elias;
  elias.add_offset(random_array, len, 1);
  auto mid = high_resolution_clock::now();
  elias.add_offset(random_array, len, -1);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  auto mid_duration = duration_cast<microseconds>(mid - start);
  auto decomp_duration = duration_cast<microseconds>(stop - mid);
  std::cout << "Total: " << duration.count() << " forward: " << mid_duration.count()
            << " backward: " << decomp_duration.count() << std::endl;
  delete[] random_array;
}

TEST(Helpers_add_offset, CheckValues) {
#define SIZE 9
  std::size_t size = SIZE;
  long input[SIZE] = {
      1, -1, -5, 5, -100, 100, 10000, -10000, -5,
  };
  long input_gt[SIZE] = {1, -1, -5, 5, -100, 100, 10000, -10000, -5};
  long results_gt[SIZE] = {2, 0, -4, 6, -99, 101, 10001, -9999, -4};
  compc::EliasOmega<long> elias;
  elias.add_offset(input, size, 1);
  for (int i = 0; i < SIZE; i++) {
    ASSERT_EQ(input[i], results_gt[i]); // comparing values
  }
  elias.add_offset(input, size, -1);
  for (int i = 0; i < SIZE; i++) {
    // std::cout << input[i] << std::endl;
    ASSERT_EQ(input[i], input_gt[i]); // comparing values
  }
}

TEST(Elias_Omega_DecompCompEQTestNegativeLong, CheckValues) {
  std::size_t size = 10;
  long input[10] = {1, 3, 2000, 2, 50, 1, 25345, 11, 10000, 1};
  compc::EliasOmega<long> elias{1, true};
  ASSERT_EQ(elias.offset, 1);
  uint8_t* comp = elias.compress(input, size);
  long* output = elias.decompress(comp, size, 10);
  for (int i = 0; i < 10; i++) {
    ASSERT_EQ(output[i], input[i]); // comparing values
  }
  delete[] comp;
  delete[] output;
}

// TODO: For offset and mapping to numbers we are not doing an overflow check.
// The above test fails for short.

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
