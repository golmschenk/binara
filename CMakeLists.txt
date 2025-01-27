cmake_minimum_required(VERSION 3.28)
project(binara)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

set(BINARA_EXE_AND_EXT_SHARED_SOURCES
        src/mcmc_wrapper.c
        src/likelihood.c
        src/likelihood.h
        src/util.c
        src/util.h
        src/python_interrupt_handling.cpp
        src/python_interrupt_handling.h
        src/random_generator.cpp
        src/random_generator.h
)

add_executable(binara_exe
        ${BINARA_EXE_AND_EXT_SHARED_SOURCES}
)
target_link_libraries(binara_exe PRIVATE OpenMP::OpenMP_CXX)

# When installed as a Python package (via pip using scikit-build-core).
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

    target_link_libraries(binara_ext PRIVATE OpenMP::OpenMP_CXX)
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
