if(CHARGEFW_PYTHON_EXECUTABLE)
    set(Python_EXECUTABLE "${CHARGEFW_PYTHON_EXECUTABLE}")
endif()
find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
set(CHARGEFW_RESOLVED_PYTHON_EXECUTABLE "${Python_EXECUTABLE}" PARENT_SCOPE)
if(SKBUILD)
    set(CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY "chargefw")
else()
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c
                "import sysconfig; print(sysconfig.get_path('platlib', vars=dict(base='', platbase='')).lstrip('/\\\\'))"
        OUTPUT_VARIABLE CHARGEFW_PYTHON_SITEARCH_RELATIVE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )
    if(CHARGEFW_PYTHON_SITEARCH_RELATIVE STREQUAL "")
        message(FATAL_ERROR "Could not determine the Python platform package installation directory")
    endif()
    set(CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY
        "${CHARGEFW_PYTHON_SITEARCH_RELATIVE}/chargefw")
endif()
set(CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY
    "${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}" PARENT_SCOPE)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_VARIABLE nanobind_ROOT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)
list(PREPEND CMAKE_PREFIX_PATH "${nanobind_ROOT}")
find_package(nanobind CONFIG REQUIRED)

set(CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY "${PROJECT_BINARY_DIR}/python/chargefw")
file(MAKE_DIRECTORY "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}")
file(REMOVE_RECURSE "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/adapters")
set(CHARGEFW_PYTHON_IO_DIRECTORY
    "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/io")
file(MAKE_DIRECTORY "${CHARGEFW_PYTHON_IO_DIRECTORY}")
set(CHARGEFW_PYTHON_STUB_DIRECTORY "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/_chargefw")
set(CHARGEFW_PYTHON_PARAMETER_DIRECTORY
    "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/_data/parameters")
file(REMOVE_RECURSE
    "${CHARGEFW_PYTHON_STUB_DIRECTORY}"
    "${CHARGEFW_PYTHON_PARAMETER_DIRECTORY}"
)
file(MAKE_DIRECTORY "${CHARGEFW_PYTHON_STUB_DIRECTORY}"
                    "${CHARGEFW_PYTHON_PARAMETER_DIRECTORY}")

file(GLOB CHARGEFW_PYTHON_PARAMETER_SOURCES CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/data/parameters/*.json")
foreach(parameter_source IN LISTS CHARGEFW_PYTHON_PARAMETER_SOURCES)
    get_filename_component(parameter_name "${parameter_source}" NAME)
    configure_file("${parameter_source}"
                   "${CHARGEFW_PYTHON_PARAMETER_DIRECTORY}/${parameter_name}" COPYONLY)
endforeach()

set(CHARGEFW_PYTHON_STUB_SOURCES
    __init__.pyi
    adapters.pyi
    calculation.pyi
    core.pyi
    methods.pyi
    parameters.pyi
)
foreach(stub_source IN LISTS CHARGEFW_PYTHON_STUB_SOURCES)
    configure_file("chargefw/_chargefw/${stub_source}"
                   "${CHARGEFW_PYTHON_STUB_DIRECTORY}/${stub_source}" COPYONLY)
endforeach()

set(CHARGEFW_PYTHON_PACKAGE_SOURCES
    __init__.py
    _arrays.py
    _catalog.py
    _elements.py
    _calculation_options.py
    _calculation_values.py
    _methods.py
    _parameters.py
    _payloads.py
    _resources.py
    _types.py
    calculation.py
    charges.py
    core.py
)
foreach(package_source IN LISTS CHARGEFW_PYTHON_PACKAGE_SOURCES)
    configure_file("chargefw/${package_source}"
                   "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/${package_source}" COPYONLY)
endforeach()
set(CHARGEFW_PYTHON_IO_SOURCES
    __init__.py
    gemmi.py
    rdkit.py
)
foreach(io_source IN LISTS CHARGEFW_PYTHON_IO_SOURCES)
    configure_file("chargefw/io/${io_source}"
                   "${CHARGEFW_PYTHON_IO_DIRECTORY}/${io_source}" COPYONLY)
endforeach()
configure_file(chargefw/py.typed "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/py.typed" COPYONLY)

nanobind_add_module(chargefw_python NB_STATIC NB_SUPPRESS_WARNINGS
    src/adapters.cpp
    src/calculation.cpp
    src/core.cpp
    src/methods.cpp
    src/module.cpp
    src/parameters.cpp
)
set_target_properties(chargefw_python PROPERTIES
    OUTPUT_NAME _chargefw
    LIBRARY_OUTPUT_DIRECTORY "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}"
    INSTALL_RPATH "${CHARGEFW_RUNTIME_ORIGIN}"
    INSTALL_RPATH_USE_LINK_PATH ON
)
# The extension contains the ChargeFW library itself, so a wheel needs no separate ChargeFW library and
# Windows needs no exported symbols from one.
target_link_libraries(chargefw_python PRIVATE chargefw_core_objects)
target_include_directories(chargefw_python PRIVATE "${PROJECT_SOURCE_DIR}/src")
chargefw_prefer_runtime_path(chargefw_python)

if(CHARGEFW_ENABLE_WARNINGS)
    chargefw_enable_warnings(chargefw_python)
endif()

if(CHARGEFW_ENABLE_ASAN)
    chargefw_enable_sanitizer(chargefw_python address)
endif()

if(CHARGEFW_ENABLE_UBSAN)
    chargefw_enable_sanitizer(chargefw_python undefined)
endif()

install(TARGETS chargefw_python
    LIBRARY DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY} COMPONENT python
    RUNTIME DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY} COMPONENT python
)
foreach(package_source IN LISTS CHARGEFW_PYTHON_PACKAGE_SOURCES)
    install(FILES "chargefw/${package_source}"
            DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY} COMPONENT python)
endforeach()
foreach(io_source IN LISTS CHARGEFW_PYTHON_IO_SOURCES)
    install(FILES "chargefw/io/${io_source}"
            DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}/io
            COMPONENT python)
endforeach()
install(FILES chargefw/py.typed
        DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}
        COMPONENT python)
install(DIRECTORY "${PROJECT_SOURCE_DIR}/python/chargefw/_chargefw"
        DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}
        COMPONENT python
        FILES_MATCHING PATTERN "*.pyi")
install(DIRECTORY "${PROJECT_SOURCE_DIR}/data/parameters/"
        DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}/_data/parameters
        COMPONENT python
        FILES_MATCHING PATTERN "*.json")
# Fetched dependency libraries are installed beside the extension, which finds them through its
# runtime search path.
if(TARGET gemmi_cpp)
    install(TARGETS gemmi_cpp
        LIBRARY DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY} COMPONENT python
        RUNTIME DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY} COMPONENT python
    )
endif()

if(TARGET tbb)
    install(TARGETS tbb
        LIBRARY DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY} COMPONENT python
        RUNTIME DESTINATION ${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY} COMPONENT python
    )
endif()
