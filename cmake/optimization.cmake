# Portable x86-64 builds, including the wheels and the CLI image, target x86-64-v3 (AVX2 and FMA,
# Intel Haswell and AMD Excavator or newer). MSVC has no host-native option, so it always uses
# x86-64-v3 on x86-64. Other architectures use the compiler default.
set(CHARGEFW_X86_64_V3 OFF)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$"
   AND (MSVC OR NOT CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS))
  set(CHARGEFW_X86_64_V3 ON)
endif()

function(chargefw_enable_optimizations target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    get_target_property(target_type ${target} TYPE)

    if(target_type MATCHES "^(SHARED|OBJECT)_LIBRARY$" AND UNIX AND NOT APPLE)
      # ChargeFW does not support runtime replacement of its exported
      # functions. Allow internal calls in the shared library and the Python
      # extension to be optimized like a static build; without this, semantic
      # interposition costs about 40% here. Only ELF platforms interpose
      # symbols this way.
      target_compile_options(${target}
        PRIVATE
          -fno-semantic-interposition
      )
    endif()

    if(CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS)
      target_compile_options(${target} PRIVATE -march=native)
    elseif(CHARGEFW_X86_64_V3)
      target_compile_options(${target} PRIVATE -march=x86-64-v3)
    endif()
  elseif(MSVC)
    if(CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS)
      message(WARNING "Native ChargeFW optimizations are not configured for MSVC; using x86-64-v3")
    endif()
    if(CHARGEFW_X86_64_V3)
      target_compile_options(${target} PRIVATE /arch:AVX2)
    endif()
  elseif(CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS)
    message(WARNING "Native ChargeFW optimizations are not configured for ${CMAKE_CXX_COMPILER_ID}")
  endif()
endfunction()
