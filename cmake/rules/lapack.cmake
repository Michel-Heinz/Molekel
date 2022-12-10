### sets linalg as a list of linear algebra libraries

# adds Amolqc/cmake dir to path that is searched for FindPackage.cmake files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# includes FindPackageHandleStandardArgs.cmake for custom FindPackage.cmake files
include(FindPackageHandleStandardArgs)

if(NOT DEFINED ENV{MATHLIBS})
    # using FindLAPACK.cmake file provided by cmake
    find_package(LAPACK REQUIRED)
    set(linalg
            ${LAPACK_LIBRARIES})
else()
    message("found MATHLIBS $ENV{MATHLIBS}")
    # using self-written FindMATHLIBS.cmake file
    find_package(MATHLIBS REQUIRED)
    set(linalg
            ${BLAS_LIBRARY}
            ${LAPACK_LIBRARY}
            )
endif()
