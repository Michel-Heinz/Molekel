cmake_minimum_required(VERSION 3.10)

add_compile_options(-ffree-line-length-1000)

add_fortran_library(utils STATIC
        parsing.f90
        libboys_data.f90
        error.f90
        string_utility.f90
        machine.f90
        kinds.f90
        math.f90
        global.f90
        elements.f90
        )

target_link_libraries(utils
        src)