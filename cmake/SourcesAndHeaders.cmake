set(sources src/elias_gamma.cpp src/elias_delta.cpp src/elias_omega.cpp)

set(exe_sources src/main.cpp ${sources})

set(headers
    include/compintc/compressor.hpp include/compintc/elias_base.hpp
    include/compintc/elias_gamma.hpp include/compintc/elias_delta.hpp
    include/compintc/elias_omega.hpp include/compintc/helpers.hpp)

set(test_sources src/elias_gamma_test.cpp src/elias_delta_test.cpp
                 src/elias_omega_test.cpp)
