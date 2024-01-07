#
# Print a message only if the `VERBOSE_OUTPUT` option is on
#

function(verbose_message content)
  if(${PROJECT_NAME}_VERBOSE_OUTPUT)
    message(STATUS ${content})
  endif()
endfunction()

#
# Add a target for formating the project using `clang-format` (i.e: cmake
# --build build --target clang-format)
#

function(add_clang_format_target)
  if(NOT ${PROJECT_NAME}_CLANG_FORMAT_BINARY)
    find_program(${PROJECT_NAME}_CLANG_FORMAT_BINARY clang-format)
  endif()

  if(${PROJECT_NAME}_CLANG_FORMAT_BINARY)
    if(${PROJECT_NAME}_BUILD_EXECUTABLE)
      add_custom_target(
        clang-format
        COMMAND ${${PROJECT_NAME}_CLANG_FORMAT_BINARY} -i ${exe_sources}
                ${headers}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
    elseif(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
      add_custom_target(
        clang-format
        COMMAND ${${PROJECT_NAME}_CLANG_FORMAT_BINARY} -i ${headers}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
    else()
      add_custom_target(
        clang-format
        COMMAND ${${PROJECT_NAME}_CLANG_FORMAT_BINARY} -i ${sources} ${headers}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
    endif()

    message(
      STATUS
        "Format the project using the `clang-format` target (i.e: cmake --build build --target clang-format).\n"
    )
  endif()
endfunction()

#
# Add a target for linting the project using `clang-tidy` (i.e: cmake --build
# build --target clang-tidy)
#

function(add_clang_tidy_target)
  if(NOT ${PROJECT_NAME}_CLANG_TIDY_BINARY)
    find_program(${PROJECT_NAME}_CLANG_TIDY_BINARY clang-tidy)
  endif()

  set(IGNORED_CHECKS
      -cppcoreguidelines-avoid-magic-numbers,
      -readability-implicit-bool-conversion,
      -cppcoreguidelines-pro-bounds-pointer-arithmetic,
      -readability-magic-numbers,
      -altera-struct-pack-align,
      -misc-non-private-member-variables-in-classes,
      -cppcoreguidelines-init-variables, # Broken at the moment
      -altera-id-dependent-backward-branch,
      -altera-unroll-loops, # currently not supported in gcc
      -bugprone-easily-swappable-parameters,
      -bugprone-branch-clone,
      -readability-function-cognitive-complexity)
  message(IGNORED_CHECKS=${IGNORED_CHECKS})

  if(${PROJECT_NAME}_CLANG_TIDY_BINARY)
    if(${PROJECT_NAME}_BUILD_EXECUTABLE)
      add_custom_target(
        clang-tidy
        COMMAND
          ${${PROJECT_NAME}_CLANG_TIDY_BINARY} -p
          ${CMAKE_BINARY_DIR}/compile_commands.json -checks="${IGNORED_CHECKS}"
          ${exe_sources} ${headers}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
    elseif(${PROJECT_NAME}_BUILD_HEADERS_ONLY)
      add_custom_target(
        clang-tidy
        COMMAND
          ${${PROJECT_NAME}_CLANG_TIDY_BINARY} -p
          ${CMAKE_BINARY_DIR}/compile_commands.json -checks="${IGNORED_CHECKS}"
          ${headers}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
    else()
      add_custom_target(
        clang-tidy
        COMMAND
          ${${PROJECT_NAME}_CLANG_TIDY_BINARY} -p
          ${CMAKE_BINARY_DIR}/compile_commands.json -checks="${IGNORED_CHECKS}"
          ${sources} ${headers}
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
    endif()

    message(
      STATUS
        "Format the project using the `clang-tidy` target (i.e: cmake --build build --target clang-tidy).\n"
    )
  endif()
endfunction()
