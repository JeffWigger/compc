cmake_minimum_required(VERSION 3.15)

#
# Project details
#

# set(CMAKE_C_COMPILER "clang")

# set(CMAKE_CXX_COMPILER "clang++") if openmp does not work! set
# (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

project(
  "compintc"
  VERSION 0.1.0
  LANGUAGES CXX)

#
# Set project options
#

include(cmake/StandardSettings.cmake)
include(cmake/StaticAnalyzers.cmake)
include(cmake/Utils.cmake)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Started CMake for ${PROJECT_NAME} v${PROJECT_VERSION}...\n")

if(UNIX)
  add_compile_options(
    "$<$<CONFIG:DEBUG>:-D_DEBUG>") # this will allow to use same _DEBUG macro
                                   # available in both Linux as well as Windows
                                   # - MSCV environment. Easy to put Debug
                                   # specific code.
endif(UNIX)

# User defined cache variables

set(SANATIZE_FLAG
    "address"
    CACHE STRING "Option for -fsanatize")

#
# Setup alternative names
#

if(${PROJECT_NAME}_USE_ALT_NAMES)
  string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)
  string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)
else()
  set(PROJECT_NAME_LOWERCASE ${PROJECT_NAME})
  set(PROJECT_NAME_UPPERCASE ${PROJECT_NAME})
endif()

#
# Prevent building in the source directory
#

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

#
# Enable package managers
#

include(cmake/Conan.cmake)
include(cmake/Vcpkg.cmake)

#
# Create library, setup header and source files
#

# Find all headers and implementation files
include(cmake/SourcesAndHeaders.cmake)

if(${PROJECT_NAME}_BUILD_EXECUTABLE)
  add_executable(${PROJECT_NAME} ${exe_sources})

  if(${PROJECT_NAME}_VERBOSE_OUTPUT)
    verbose_message("Found the following sources:")
    foreach(source IN LISTS exe_sources)
      verbose_message("* ${source}")
    endforeach()
  endif()

  if(${PROJECT_NAME}_ENABLE_UNIT_TESTING)
    add_library(${PROJECT_NAME}_LIB ${headers} ${sources})

    if(${PROJECT_NAME}_VERBOSE_OUTPUT)
      verbose_message("Found the following headers:")
      foreach(header IN LISTS headers)
        verbose_message("* ${header}")
      endforeach()
    endif()
  endif()
elseif(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  add_library(${PROJECT_NAME} INTERFACE)

  if(${PROJECT_NAME}_VERBOSE_OUTPUT)
    verbose_message("Found the following headers:")
    foreach(header IN LIST headers)
      verbose_message("* ${header}")
    endforeach()
  endif()
else()
  add_library(${PROJECT_NAME} ${headers} ${sources})

  if(${PROJECT_NAME}_VERBOSE_OUTPUT)
    verbose_message("Found the following sources:")
    foreach(source IN LISTS sources)
      verbose_message("* ${source}")
    endforeach()
    verbose_message("Found the following headers:")
    foreach(header IN LISTS headers)
      verbose_message("* ${header}")
    endforeach()
  endif()
endif()

set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES ARCHIVE_OUTPUT_DIRECTORY
             "${CMAKE_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE}"
             LIBRARY_OUTPUT_DIRECTORY
             "${CMAKE_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE}"
             RUNTIME_OUTPUT_DIRECTORY
             "${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE}")
if(${PROJECT_NAME}_BUILD_EXECUTABLE AND ${PROJECT_NAME}_ENABLE_UNIT_TESTING)
  set_target_properties(
    ${PROJECT_NAME}_LIB
    PROPERTIES ARCHIVE_OUTPUT_DIRECTORY
               "${CMAKE_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE}"
               LIBRARY_OUTPUT_DIRECTORY
               "${CMAKE_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE}"
               OUTPUT_NAME ${PROJECT_NAME})
endif()

message(STATUS "Added all header and implementation files.\n")

#
# Set the project build flags
#

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  message(STATUS "Using -O3 Flag\n")
  # May need  -Xclang for clang could conflict with flags set by tests. For
  # clang needs libomp-dev dependency
  target_compile_options(${PROJECT_NAME} PRIVATE -O3)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Sanatize")
  message(STATUS "Sanatize + Using -O3 Flag\n")
  target_compile_options(${PROJECT_NAME} PRIVATE -O3 -g) # without -g we do not
                                                         # get the line numbers.
  target_compile_options(${PROJECT_NAME} PUBLIC -fsanitize=${SANATIZE_FLAG}
  )# -fsanitize-memory-track-origins -fsanitize=undefined -fsanitize=memory
  # is not working
  target_link_options(${PROJECT_NAME} PUBLIC -fsanitize=${SANATIZE_FLAG}
  )# -fsanitize-memory-track-origins -fsanitize=undefined -fsanitize=memory
  if("${SANATIZE_FLAG}" STREQUAL "memory")
    target_compile_options(${PROJECT_NAME}
                           PUBLIC -fsanitize-memory-track-origins=2)
    target_link_options(${PROJECT_NAME} PUBLIC
                        -fsanitize-memory-track-origins=2)
  endif()
  # add_link_options(-fsanitize=${SANATIZE_FLAG})
else()
  message(STATUS "Using -O0 Flag\n")
  target_compile_options(${PROJECT_NAME} PRIVATE -O0 -g)
endif()

#
# Set the project standard and warnings
#

if(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
else()
  target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

  if(${PROJECT_NAME}_BUILD_EXECUTABLE AND ${PROJECT_NAME}_ENABLE_UNIT_TESTING)
    target_compile_features(${PROJECT_NAME}_LIB PUBLIC cxx_std_17)
  endif()
endif()
include(cmake/CompilerWarnings.cmake)
set_project_warnings(${PROJECT_NAME})

verbose_message(
  "Applied compiler warnings. Using standard ${CMAKE_CXX_STANDARD}.\n")

#
# Enable Doxygen
#

include(cmake/Doxygen.cmake)

#
# Model project dependencies
#
# May need  -Xclang for clang could conflict with flags set by tests. For clang
# needs libomp-dev dependency
find_package(OpenMP REQUIRED)
target_link_libraries(
  ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX) # Needs to be public otherwise it
                                             # does not work sometimes
# Identify and link with the specific "packages" the project uses
# find_package(package_name package_version REQUIRED package_type
# [other_options]) target_link_libraries( ${PROJECT_NAME} PUBLIC dependency1 ...
# PRIVATE dependency2 ... ${PROJECT_NAME}_PROJECT_OPTIONS
# ${PROJECT_NAME}_PROJECT_WARNINGS ) if(${PROJECT_NAME}_BUILD_EXECUTABLE AND
# ${PROJECT_NAME}_ENABLE_UNIT_TESTING) target_link_libraries(
# ${PROJECT_NAME}_LIB PUBLIC dependency1 ... ) endif()

# For Windows, it is necessary to link with the MultiThreaded library. Depending
# on how the rest of the project's dependencies are linked, it might be
# necessary to change the line to statically link with the library.
#
# This is done as follows:
#
# set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
#
# On Linux and Mac this variable is ignored. If any issues rise from it, try
# commenting it out and letting CMake decide how to link with it.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")

verbose_message("Successfully added all dependencies and linked against them.")

#
# Set the build/user include directories
#

# Allow usage of header files in the `src` directory, but only for utilities
if(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
  target_include_directories(
    ${PROJECT_NAME}
    INTERFACE $<INSTALL_INTERFACE:include>
              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
else()
  target_include_directories(
    ${PROJECT_NAME}
    PUBLIC $<INSTALL_INTERFACE:include>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  if(${PROJECT_NAME}_BUILD_EXECUTABLE AND ${PROJECT_NAME}_ENABLE_UNIT_TESTING)
    target_include_directories(
      ${PROJECT_NAME}_LIB
      PUBLIC $<INSTALL_INTERFACE:include>
             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  endif()
endif()

message(STATUS "Finished setting up include directories.")

#
# Provide alias to library for
#

if(${PROJECT_NAME}_BUILD_EXECUTABLE)
  add_executable(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
else()
  add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
endif()

verbose_message("Project is now aliased as ${PROJECT_NAME}::${PROJECT_NAME}.\n")

#
# Format the project using the `clang-format` target (i.e: cmake --build build
# --target clang-format)
#

add_clang_format_target()

#
# Format the project using the `clang-format` target (i.e: cmake --build build
# --target clang-tidy)
#

add_clang_tidy_target()

#
# Install library for easy downstream inclusion
#

include(GNUInstallDirs)
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES
  DESTINATION include
  PUBLIC_HEADER DESTINATION include)

install(
  EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

#
# Add version header
#

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/version.hpp.in
               include/${PROJECT_NAME_LOWERCASE}/version.hpp @ONLY)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME_LOWERCASE}/version.hpp
  DESTINATION include/${PROJECT_NAME_LOWERCASE})

#
# Install the `include` directory
#

install(DIRECTORY include/${PROJECT_NAME_LOWERCASE} DESTINATION include)

verbose_message(
  "Install targets successfully built. Install with `cmake --build <build_directory> --target install --config <build_config>`."
)

#
# Quick `ConfigVersion.cmake` creation
#

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

#
# Generate export header if specified
#

if(${PROJECT_NAME}_GENERATE_EXPORT_HEADER)
  include(GenerateExportHeader)
  generate_export_header(${PROJECT_NAME})
  install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME_LOWERCASE}_export.h
          DESTINATION include)

  message(
    STATUS
      "Generated the export header `${PROJECT_NAME_LOWERCASE}_export.h` and installed it."
  )
endif()

message(STATUS "Finished building requirements for installing the package.\n")

#
# Unit testing setup
#

if(${PROJECT_NAME}_ENABLE_UNIT_TESTING)
  enable_testing()
  message(
    STATUS
      "Build unit tests for the project. Tests should always be found in the test folder\n"
  )
  add_subdirectory(test)
endif()

get_target_property(MAIN_CFLAGS ${PROJECT_NAME} COMPILE_OPTIONS)
# also see: COMPILE_DEFINITIONS INCLUDE_DIRECTORIES
message("-- Target compiler flags are: ${MAIN_CFLAGS}")
