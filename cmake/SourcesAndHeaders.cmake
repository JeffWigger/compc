set(sources
    src/elias_gamma.cpp
)

set(exe_sources
		src/main.cpp 
		${sources}
)

set(headers
    include/compc/compressor.hpp
    include/compc/elias_base.hpp
    include/compc/elias_gamma.hpp
    include/compc/helpers.hpp 
)

set(test_sources
  src/elias_test.cpp
)
