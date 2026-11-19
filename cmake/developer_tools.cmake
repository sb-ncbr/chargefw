function(chargefw_setup_developer_tools)
  if(CHARGEFW_ENABLE_CCACHE)
    find_program(CCACHE_PROGRAM ccache)

    if(CCACHE_PROGRAM)
      message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
      set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}" PARENT_SCOPE)
    else()
      message(STATUS "ccache not found")
    endif()
  endif()

  if(CHARGEFW_ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE clang-tidy REQUIRED)

    # Eigen's temporary buffers come from alloca or malloc, and the analyzer pairs one branch's
    # allocation with the other's release. Analyze the heap-only form; the build is unchanged.
    set(
      CMAKE_CXX_CLANG_TIDY
      "${CLANG_TIDY_EXE}"
      "--header-filter=^${PROJECT_SOURCE_DIR}/(include|src|apps|tests)/"
      "--exclude-header-filter=^(${PROJECT_BINARY_DIR}|${PROJECT_SOURCE_DIR}/build)/"
      "--extra-arg=-DEIGEN_NO_ALLOCA"
      PARENT_SCOPE
    )
  else()
    set(CMAKE_CXX_CLANG_TIDY "" PARENT_SCOPE)
  endif()
endfunction()
