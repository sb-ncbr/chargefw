# Installed dependencies outside the system library paths are referenced by their absolute
# directories; libraries built within this project are found relative to the installation.
set_target_properties(chargefw_core
        PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        INSTALL_RPATH_USE_LINK_PATH ON
)

# Prefer the installed dependency copies over same-named libraries from
# LD_LIBRARY_PATH. This is important for fetched dependencies such as Gemmi,
# whose unversioned SONAME can collide with another installation.
target_link_options(chargefw_core PRIVATE "LINKER:--disable-new-dtags")

if(TARGET chargefw_cli)
    set_target_properties(chargefw_cli
            PROPERTIES
            INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
            INSTALL_RPATH_USE_LINK_PATH ON
    )
    target_link_options(chargefw_cli PRIVATE "LINKER:--disable-new-dtags")
endif()

install(
        TARGETS chargefw_core
        EXPORT chargefwTargets
        FILE_SET HEADERS
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        ARCHIVE
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME
        DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
        DIRECTORY
        ${PROJECT_SOURCE_DIR}/data/parameters
        DESTINATION
        ${CMAKE_INSTALL_DATADIR}/chargefw
        FILES_MATCHING
        PATTERN "*.json"
)

# Keep the public guide hierarchy intact so installed Markdown links remain valid.
install(
        FILES
        ${PROJECT_SOURCE_DIR}/README.md
        ${PROJECT_SOURCE_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DOCDIR}/chargefw
        COMPONENT documentation
)
install(
        DIRECTORY
        ${PROJECT_SOURCE_DIR}/licenses
        DESTINATION ${CMAKE_INSTALL_DOCDIR}/chargefw
        COMPONENT documentation
)
install(
        DIRECTORY ${PROJECT_SOURCE_DIR}/docs/
        DESTINATION ${CMAKE_INSTALL_DOCDIR}/chargefw/docs
        COMPONENT documentation
        FILES_MATCHING
        PATTERN "*.md"
        PATTERN "*.py"
        PATTERN "*.sdf"
        PATTERN "*.cif"
        PATTERN "conf.py" EXCLUDE
        PATTERN "_ext" EXCLUDE
        PATTERN "_generated" EXCLUDE
        PATTERN "parameter-sets.md" EXCLUDE
)
# The documentation site generates the parameter-set tables; installations get a page without them.
install(
        FILES ${PROJECT_SOURCE_DIR}/cmake/installed-parameter-sets.md
        DESTINATION ${CMAKE_INSTALL_DOCDIR}/chargefw/docs/reference
        RENAME parameter-sets.md
        COMPONENT documentation
)
install(
        FILES
        ${PROJECT_SOURCE_DIR}/schemas/molecule-input-1.0.schema.json
        ${PROJECT_SOURCE_DIR}/schemas/result-1.0.schema.json
        DESTINATION ${CMAKE_INSTALL_DOCDIR}/chargefw/schemas
        COMPONENT documentation
)

install(
        FILES
        ${PROJECT_BINARY_DIR}/generated/chargefw/config.h
        DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}/chargefw
)

configure_package_config_file(
        ${PROJECT_SOURCE_DIR}/cmake/chargefwConfig.cmake.in
        ${PROJECT_BINARY_DIR}/chargefwConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chargefw
)
write_basic_package_version_file(
        ${PROJECT_BINARY_DIR}/chargefwConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)
install(
        EXPORT chargefwTargets
        FILE chargefwTargets.cmake
        NAMESPACE chargefw::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chargefw
)
install(
        FILES
        ${PROJECT_BINARY_DIR}/chargefwConfig.cmake
        ${PROJECT_BINARY_DIR}/chargefwConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chargefw
)

if(CHARGEFW_TBB_FETCHED AND TARGET tbb)
    set_target_properties(tbb PROPERTIES EXPORT_NAME _tbb_runtime)
    install(TARGETS tbb
            EXPORT chargefwTargets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()

if(CHARGEFW_GEMMI_FETCHED AND TARGET gemmi_cpp)
    install(TARGETS gemmi_cpp
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()
