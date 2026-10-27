function(chargefw_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4)
    return()
  endif()

  target_compile_options(${target}
    PRIVATE
      -Wall
      -Wextra
      -Wpedantic
      -Wconversion
      -Wsign-conversion
      -Wshadow
      -Wnon-virtual-dtor
      -Wold-style-cast
      -Woverloaded-virtual
      -Wnull-dereference
      -Wdouble-promotion
      -Wformat=2
  )

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC emits false-positive null-dereference warnings from Eigen's LU internals in optimized
    # builds; keep the warning enabled for other compilers and project code where it is useful.
    target_compile_options(${target} PRIVATE -Wno-null-dereference)
  endif()
endfunction()
