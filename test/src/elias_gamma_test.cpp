#include "compintc/elias_gamma.hpp"
#include "helpers.hpp"
#include <bitset>
#include <chrono>
#include <cstdint>
#include <gtest/gtest.h>
using namespace std::chrono;

TEST(Elias_Gamma_DecompCompEQTestLong, CheckValues) {
  std::size_t size = 10;
  long input[10] = {1, 3, 2000, 2, 50, 1, 25345, 11, 10000000, 1};
  compc::EliasGamma<long> elias;
  std::cout << "size: " << sizeof(input) << std::endl;
  std::unique_ptr<uint8_t[]> comp = elias.compress(input, size);

  std::cout << "binary size: " << size << std::endl;
  for (std::size_t i = 0; i < size; i++)
    std::cout << std::bitset<8>(comp[i]) << std::endl;
  std::unique_ptr<long[]> output = elias.decompress(comp.get(), size, 10);
  // ASSERT_EQ(size, 21);
  for (std::size_t i = 0; i < 10; i++) {
    ASSERT_EQ(output[i], input[i]); // comparing values
  }
}

TEST(Elias_Gamma_SpeedTestLong, CheckValues) {
  std::size_t len = 500000;
  std::size_t len_copy = len;
  auto random_array = compc_test::get_random_array<long>(len);

  auto start = high_resolution_clock::now();
  compc::EliasGamma<long> elias;
  std::unique_ptr<uint8_t[]> comp = elias.compress(random_array.get(), len);
  std::cout << "binary size: " << len << std::endl;
  auto mid = high_resolution_clock::now();
  std::unique_ptr<long[]> output = elias.decompress(comp.get(), len, len_copy);
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
}

TEST(Elias_Gamma_DecompCompEQTestInt, CheckValues) {
  std::size_t size = 10;
  int input[10] = {1, 3, 2000, 2, 50, 1, 25345, 11, 1000000, 1};
  compc::EliasGamma<int> elias;
  std::cout << "size: " << sizeof(input) << std::endl;
  std::unique_ptr<uint8_t[]> comp = elias.compress(input, size);

  std::unique_ptr<int[]> output = elias.decompress(comp.get(), size, 10);
  // ASSERT_EQ(size, 21);
  for (std::size_t i = 0; i < 10; i++) {
    ASSERT_EQ(output[i], input[i]); // comparing values
  }
}

TEST(Elias_Gamma_DecompCompEQTestShort, CheckValues) {
  std::size_t size = 10;
  short input[10] = {1, 3, 2000, 2, 50, 1, 25345, 11, 10000, 1};
  compc::EliasGamma<short> elias;
  std::cout << "size: " << sizeof(input) << std::endl;
  std::unique_ptr<uint8_t[]> comp = elias.compress(input, size);

  std::unique_ptr<short[]> output = elias.decompress(comp.get(), size, 10);
  // ASSERT_EQ(size, 21);
  for (std::size_t i = 0; i < 10; i++) {
    ASSERT_EQ(output[i], input[i]); // comparing values
  }
}

TEST(Helpers_transform_to_natural_numbers, CheckValues) {
  std::size_t size = 8;
  long input[8] = {1, -1, -5, 5, -100, 100, 10000, -10000};
  long input_gt[8] = {1, -1, -5, 5, -100, 100, 10000, -10000};
  long results_gt[8] = {2, 1, 9, 10, 199, 200, 20000, 19999};
  compc::EliasGamma<long> elias;
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
  compc::EliasGamma<long> elias;
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
  compc::EliasGamma<long> elias;
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
  auto random_array = compc_test::get_random_array<long>(len);

  auto start = high_resolution_clock::now();
  compc::EliasGamma<long> elias;
  elias.transform_to_natural_numbers(random_array.get(), len);
  auto mid = high_resolution_clock::now();
  elias.transform_to_natural_numbers_reverse(random_array.get(), len);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  auto mid_duration = duration_cast<microseconds>(mid - start);
  auto decomp_duration = duration_cast<microseconds>(stop - mid);
  std::cout << "Total: " << duration.count() << " forward: " << mid_duration.count()
            << " backward: " << decomp_duration.count() << std::endl;
}

TEST(Helpers_add_offset_SpeedTestLong, CheckValues) {
  std::size_t len = 500000;
  auto random_array = compc_test::get_random_array<long>(len);

  auto start = high_resolution_clock::now();
  compc::EliasGamma<long> elias;
  elias.add_offset(random_array.get(), len, 1);
  auto mid = high_resolution_clock::now();
  elias.add_offset(random_array.get(), len, -1);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  auto mid_duration = duration_cast<microseconds>(mid - start);
  auto decomp_duration = duration_cast<microseconds>(stop - mid);
  std::cout << "Total: " << duration.count() << " forward: " << mid_duration.count()
            << " backward: " << decomp_duration.count() << std::endl;
}

TEST(Helpers_add_offset, CheckValues) {
#define SIZE 9
  std::size_t size = SIZE;
  long input[SIZE] = {
      1, -1, -5, 5, -100, 100, 10000, -10000, -5,
  };
  long input_gt[SIZE] = {1, -1, -5, 5, -100, 100, 10000, -10000, -5};
  long results_gt[SIZE] = {2, 0, -4, 6, -99, 101, 10001, -9999, -4};
  compc::EliasGamma<long> elias;
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

TEST(Elias_Gamma_DecompCompEQTestNegativeLong, CheckValues) {
  std::size_t size = 10;
  long input[10] = {1, 3, 2000, 2, 50, 1, 25345, 11, 10000, 1};
  compc::EliasGamma<long> elias{1, true};
  ASSERT_EQ(elias.offset, 1);
  std::unique_ptr<uint8_t[]> comp = elias.compress(input, size);
  std::unique_ptr<long[]> output = elias.decompress(comp.get(), size, 10);
  for (std::size_t i = 0; i < 10; i++) {
    ASSERT_EQ(output[i], input[i]); // comparing values
  }
}

// TODO: For offset and mapping to numbers we are not doing an overflow check.
// The above test fails for short.

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
