cmake_minimum_required(VERSION 3.10)

add_fortran_library(eval STATIC
        RHF.f90
        initEval.f90
        dipolMoment.f90
        )
add_dependencies (eval integrals)
target_link_libraries(eval
        utils
        ${linalg}
        integrals
        src)