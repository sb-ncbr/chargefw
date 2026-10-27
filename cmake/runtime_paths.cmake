# Installed binaries find ChargeFW's own libraries relative to themselves. Windows has no runtime
# search paths; it loads DLLs from the directory of the executable or extension module instead.
if(APPLE)
  set(CHARGEFW_RUNTIME_ORIGIN "@loader_path")
else()
  set(CHARGEFW_RUNTIME_ORIGIN "$ORIGIN")
endif()

# Prefer the bundled dependency copies over same-named libraries from LD_LIBRARY_PATH. This is
# important for fetched dependencies such as Gemmi, whose unversioned SONAME can collide with another
# installation. Only ELF platforms distinguish RPATH from RUNPATH.
function(chargefw_prefer_runtime_path target)
  if(UNIX AND NOT APPLE)
    target_link_options(${target} PRIVATE "LINKER:--disable-new-dtags")
  endif()
endfunction()
