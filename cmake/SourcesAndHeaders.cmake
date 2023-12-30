set(sources
    src/elias_gamma.cpp
    src/elias_delta.cpp
    src/elias_omega.cpp
)

set(exe_sources
		src/main.cpp 
		${sources}
)

set(headers
    include/compc/compressor.hpp
    include/compc/elias_base.hpp
    include/compc/elias_gamma.hpp
    include/compc/elias_delta.hpp
    include/compc/elias_omega.hpp
    include/compc/helpers.hpp 
)

set(test_sources
    src/elias_gamma_test.cpp
    src/elias_delta_test.cpp
    src/elias_omega_test.cpp
)
