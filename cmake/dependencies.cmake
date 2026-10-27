include(FetchContent)

function(chargefw_setup_dependencies)
    set(CHARGEFW_TBB_FETCHED OFF PARENT_SCOPE)
    set(CHARGEFW_GEMMI_FETCHED OFF PARENT_SCOPE)

    if(CHARGEFW_BUILD_CLI OR CHARGEFW_BUILD_PYTHON)
        if(NOT CHARGEFW_FETCH_DEPENDENCIES)
            find_package(CLI11 2.4 CONFIG REQUIRED)
        else()
            FetchContent_Declare(
                    cli11
                    SYSTEM
                    EXCLUDE_FROM_ALL
                    URL https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.7.2.tar.gz
                    URL_HASH SHA256=46eef3101da70852ec7af026e09d485ccee81813331c8c6052d39344443b83da
            )
            FetchContent_MakeAvailable(cli11)
        endif()
    endif()

    if(NOT CHARGEFW_FETCH_DEPENDENCIES)
        find_package(nlohmann_json 3.11.3 CONFIG REQUIRED)
        set(CHARGEFW_NLOHMANN_JSON_VERSION "${nlohmann_json_VERSION}" PARENT_SCOPE)
    else()
        set(JSON_Install ON)
        set(CHARGEFW_NLOHMANN_JSON_VERSION 3.12.0 PARENT_SCOPE)
        FetchContent_Declare(
                nlohmann_json
                SYSTEM
                URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
                URL_HASH SHA256=42f6e95cad6ec532fd372391373363b62a14af6d771056dbfc86160e6dfff7aa
        )
        FetchContent_MakeAvailable(nlohmann_json)
    endif()

    if(NOT CHARGEFW_FETCH_DEPENDENCIES)
        # Eigen version files accept only one major version per request; 3.4 and 5 are supported.
        find_package(Eigen3 5 CONFIG QUIET)
        if(NOT Eigen3_FOUND)
            find_package(Eigen3 3.4 CONFIG REQUIRED)
        endif()
    else()
        FetchContent_Declare(
                eigen
                SYSTEM
                EXCLUDE_FROM_ALL
                URL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz
                URL_HASH SHA256=e9c326dc8c05cd1e044c71f30f1b2e34a6161a3b6ecf445d56b53ff1669e3dec
        )
        FetchContent_MakeAvailable(eigen)
    endif()

    if(NOT CHARGEFW_FETCH_DEPENDENCIES)
        find_package(nanoflann 1.7 CONFIG REQUIRED)
    else()
        # nanoflann currently does not set CMP0077 itself. Select modern option() behavior so these
        # embedding-project choices remain effective on the first configure.
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
        set(NANOFLANN_BUILD_EXAMPLES OFF)
        set(NANOFLANN_BUILD_TESTS OFF)
        FetchContent_Declare(
                nanoflann
                SYSTEM
                EXCLUDE_FROM_ALL
                URL https://github.com/jlblancoc/nanoflann/archive/refs/tags/2.0.1.tar.gz
                URL_HASH SHA256=bb3991ad924248d3f20953a951c34878196c11ca112d800acc341c3808443896
        )
        FetchContent_MakeAvailable(nanoflann)
        unset(CMAKE_POLICY_DEFAULT_CMP0077)
    endif()

    if(NOT CHARGEFW_FETCH_DEPENDENCIES)
        find_package(TBB 2022.1 CONFIG REQUIRED)
    else()
        # An explicit project-wide IPO choice already initializes fetched targets. Avoid oneTBB's
        # redundant all-enabled-language probe, which is unreliable when C is enabled by a sibling.
        if(DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
            set(TBB_ENABLE_IPO OFF)
        endif()
        set(BUILD_SHARED_LIBS ON)
        set(TBB_TEST OFF)
        set(TBB_STRICT OFF)
        set(TBB_INSTALL OFF)
        set(CHARGEFW_TBB_FETCHED ON PARENT_SCOPE)
        FetchContent_Declare(
                onetbb
                SYSTEM
                URL https://github.com/uxlfoundation/oneTBB/archive/refs/tags/v2023.1.0.tar.gz
                URL_HASH SHA256=191288b52e1e6b17198000b64d77d194bb65e791be46ebc606e9b091781e2070
        )
        FetchContent_MakeAvailable(onetbb)
        # A wheel stores each of a library's file names as a separate copy. Name the Python package's
        # oneTBB after its interface version only, so a single file serves the extension.
        if(NOT CHARGEFW_BUILD_LIBRARY)
            get_target_property(tbb_soversion tbb SOVERSION)
            set_target_properties(tbb PROPERTIES VERSION ${tbb_soversion})
        endif()
    endif()

    if(NOT CHARGEFW_FETCH_DEPENDENCIES)
        find_package(gemmi 0.7.4 CONFIG REQUIRED)
        set(CHARGEFW_GEMMI_VERSION "${gemmi_VERSION}" PARENT_SCOPE)
    else()
        set(CHARGEFW_GEMMI_VERSION 0.7.4 PARENT_SCOPE)
        set(BUILD_GEMMI_PROGRAM OFF)
        # Native installations ship Gemmi as a shared library, because Gemmi types appear in ChargeFW's
        # public headers. The Python package links it into the extension and the command line, so it
        # cannot meet another library of the same unversioned name in one process.
        set(BUILD_SHARED_LIBS ${CHARGEFW_BUILD_LIBRARY})
        set(INSTALL_DEV_FILES ${CHARGEFW_BUILD_LIBRARY})
        set(CHARGEFW_GEMMI_FETCHED ON PARENT_SCOPE)
        FetchContent_Declare(
                gemmi
                SYSTEM
                URL https://github.com/project-gemmi/gemmi/archive/refs/tags/v0.7.4.tar.gz
                URL_HASH SHA256=58ed9023e0f75033ebec1da630461bc6e1af661f29f0c41deed66c20315ebe2a
        )
        FetchContent_MakeAvailable(gemmi)
    endif()
endfunction()

function(chargefw_setup_test_dependencies)
    find_package(snitch 1.3.2 CONFIG QUIET)
    if(NOT TARGET snitch::snitch)
        # Snitch assertions must remain active under NDEBUG so that release
        # builds execute the same checks as debug builds.
        set(SNITCH_ENABLE ON CACHE BOOL "" FORCE)
        FetchContent_Declare(
                snitch
                SYSTEM
                EXCLUDE_FROM_ALL
                URL https://github.com/cschreib/snitch/archive/refs/tags/v1.3.2.tar.gz
                URL_HASH SHA256=5b3a0af7e5637f02ac95652f2238b89eb227b4881bf4e5ffc6935f184ff11626
        )
        FetchContent_MakeAvailable(snitch)

        if(CHARGEFW_ENABLE_CLANG_TIDY)
            # CMAKE_CXX_CLANG_TIDY also reaches FetchContent targets. Snitch is external code and
            # must not add third-party diagnostics to the project's clang-tidy output.
            set_property(TARGET snitch PROPERTY CXX_CLANG_TIDY "")
        endif()
    endif()
endfunction()
