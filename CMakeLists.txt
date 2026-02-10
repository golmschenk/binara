cmake_minimum_required(VERSION 3.28)
project(binara)

set(CMAKE_CXX_STANDARD 20)

# Find expected dependencies.
find_package(OpenMP REQUIRED)
if(POLICY CMP0167)  # Enforce using modern Boost finding method.
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost COMPONENTS REQUIRED random)

if(NOT WIN32)
    set(MATH_LIBRARY m)
else()
    set(MATH_LIBRARY "")
endif()

# Add downloadable dependencies.
include(FetchContent)
FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG        v3.4.0
)
FetchContent_MakeAvailable(tomlplusplus)

# Prepare shared source files.
set(BINARA_EXE_AND_EXT_SHARED_SOURCES
        src/mcmc_wrapper.cpp
        src/likelihood.cpp
        src/likelihood.h
        src/util.cpp
        src/util.h
        src/python_interrupt_handling.cpp
        src/python_interrupt_handling.h
        src/random_generator.cpp
        src/random_generator.h
        src/configuration.cpp
        src/configuration.h
        src/eclipsing_method_enum.cpp
        src/eclipse_method_enum.h
)

# Build executable.
add_executable(binara_exe
        ${BINARA_EXE_AND_EXT_SHARED_SOURCES}
)
target_link_libraries(binara_exe PRIVATE OpenMP::OpenMP_CXX ${MATH_LIBRARY})
target_link_libraries(binara_exe PRIVATE Boost::random tomlplusplus::tomlplusplus)

# Build Python library when installed as a Python package (via pip using scikit-build-core).
if(SKBUILD)
    find_package(Python 3.8
            REQUIRED COMPONENTS Interpreter Development.Module
            OPTIONAL_COMPONENTS Development.SABIModule)

    execute_process(
            COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
    find_package(nanobind CONFIG REQUIRED)
    nanobind_add_module(binara_ext STABLE_ABI NB_STATIC
            src/binara_ext.cpp
            ${BINARA_EXE_AND_EXT_SHARED_SOURCES}
    )

    target_link_libraries(binara_ext PRIVATE OpenMP::OpenMP_CXX ${MATH_LIBRARY})
    target_link_libraries(binara_ext PRIVATE Boost::random tomlplusplus::tomlplusplus)
    target_compile_definitions(binara_ext PRIVATE IS_PYTHON_LIBRARY_BUILD)
    nanobind_add_stub(
            binara_ext_stub
            MODULE binara_ext
            OUTPUT binara_ext.pyi
            PYTHON_PATH $<TARGET_FILE_DIR:binara_ext>
            DEPENDS binara_ext
    )

    install(TARGETS binara_ext LIBRARY DESTINATION binara)
    install(FILES "${CMAKE_BINARY_DIR}/binara_ext.pyi" DESTINATION binara)
endif()
