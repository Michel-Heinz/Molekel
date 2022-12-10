cmake_minimum_required(VERSION 3.10)

add_fortran_library(integrals STATIC
        nuclearRepulsion.f90
        integralsMain.f90
        overlapJoshuaGoings.f90
        kineticEJoshuaGoings.f90
        dipoleIntegrals.f90
        nuclearAttractionIntegralJoshuaGoings.f90
        electronRepulsionIntegralsJoshuaGoings.f90
        )
target_link_libraries(integrals
        utils
        src)