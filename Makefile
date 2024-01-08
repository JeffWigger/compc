.PHONY: install coverage test docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

# Setting default compilers
ifndef CXX
	override CXX=clang++
endif
ifndef CC
	override CC=clang
endif
$(info using ${CXX})
$(info using ${CC})



BROWSER := python -c "$$BROWSER_PYSCRIPT"
INSTALL_LOCATION := ~/.local

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

test: ## run tests quickly with ctest
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -Dcompintc_ENABLE_UNIT_TESTING=1 -DCMAKE_BUILD_TYPE="Debug" -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX) -Dcompintc_ENABLE_CODE_COVERAGE=0
	cmake --build build --config Release
	cd build/ && ctest -C Release -VV
#--gtest_filter=Elias_Gamma_DecompCompEQTestLong

test-fast: ## run tests quickly with ctest
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -Dcompintc_ENABLE_UNIT_TESTING=1 -DCMAKE_BUILD_TYPE="Release" -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX) -Dcompintc_ENABLE_CODE_COVERAGE=0
	cmake --build build --config Release
	cd build/ && ctest -C Release -VV
#--gtest_filter=Elias_Gamma_DecompCompEQTestLong

testAddress: ## run tests quickly with ctest
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -Dcompintc_ENABLE_UNIT_TESTING=1 -DCMAKE_BUILD_TYPE="Sanatize" -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX) -Dcompintc_ENABLE_CODE_COVERAGE=0 -DSANATIZE_FLAG:STRING=address
	cmake --build build --config Release
	cd build/ && ctest -C Release -VV

testMemory: ## run tests quickly with ctest
	# test not used as we deleiberately use undefined memory for performance reasons.
	# add -fsanitize-memory-track-origins for the prints to make more sense.
	rm -rf build/
	# only available with clang++
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -Dcompintc_ENABLE_UNIT_TESTING=1 -DCMAKE_BUILD_TYPE="Sanatize" -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=clang++ -Dcompintc_ENABLE_CODE_COVERAGE=0 -DSANATIZE_FLAG:STRING=memory
	cmake --build build --config Release
	cd build/ && ctest -C Release -VV

testUndefined: ## run tests quickly with ctest
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -Dcompintc_ENABLE_UNIT_TESTING=1 -DCMAKE_BUILD_TYPE="Sanatize" -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX) -Dcompintc_ENABLE_CODE_COVERAGE=0 -DSANATIZE_FLAG:STRING=undefined
	cmake --build build --config Release
	cd build/ && ctest -C Release -VV

coverage: ## check code coverage quickly GCC
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -Dcompintc_ENABLE_CODE_COVERAGE=1 -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX)
	cmake --build build --config Release
	cd build/ && ctest -C Release -VV
	cd .. && (bash -c "find . -type f -name '*.gcno' -exec gcov -pb {} +" || true)

docs: ## generate Doxygen HTML documentation, including API docs
	rm -rf docs/
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -Dcompintc_ENABLE_DOXYGEN=1 -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX)
	cmake --build build --config Release
	cmake --build build --target doxygen-docs
	$(BROWSER) docs/html/index.html

install: ## install the package to the `INSTALL_LOCATION`
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX)
	cmake --build build --config Release
	cmake --build build --target install --config Release

sanatize: ## install the package to the `INSTALL_LOCATION` with sanatization
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX) -D CMAKE_BUILD_TYPE="Sanatize"
	cmake --build build --config Release
	cmake --build build --target install --config Release

format: ## format the project sources
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION)
	cmake --build build --target clang-format

tidy: ## format the project sources
	rm -rf build/
	cmake -Bbuild -DCMAKE_INSTALL_PREFIX=$(INSTALL_LOCATION) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
	cmake --build build --target clang-tidy
