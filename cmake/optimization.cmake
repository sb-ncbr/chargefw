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
        target_compile_options(${target}
        PRIVATE
          -march=native
      )
    endif()
  elseif(CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS)
    message(WARNING "Native ChargeFW optimizations are not configured for ${CMAKE_CXX_COMPILER_ID}")
  endif()
endfunction()
