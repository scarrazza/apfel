# BASIC DEFINITIONS ========================================================================

# define minimum version of cmake
cmake_minimum_required (VERSION 3.5)

# define project name and its language
project(apfel CXX Fortran)

# define c++ standard and issue all the warning demanded by this standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -fPIC")
set(CMAKE_SHARED_LINKER_FLAGS -w)

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# VERSION ==============================================================================

set(apfel_VERSION 3.1.0)

# CONFIG SCRIPT =========================================================================

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
configure_file("${PROJECT_SOURCE_DIR}/bin/apfel-config.in" "${PROJECT_SOURCE_DIR}/bin/apfel-config")

# FINALIZE ==================================================================================

# LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config REQUIRED)
if (LHAPDF_CONFIG)
  exec_program(${LHAPDF_CONFIG}
    ARGS --cflags
    OUTPUT_VARIABLE LHAPDF_CXX_FLAGS
  )
  set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${LHAPDF_CONFIG}
    ARGS --libs
    OUTPUT_VARIABLE LHAPDF_LIBRARIES
  )
  set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARIES} CACHE STRING INTERNAL)
else(LHAPDF_CONFIG)
  message(FATAL_ERROR "LHAPDF not found!")
endif(LHAPDF_CONFIG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS}")

# include directory
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

# generate list of source files
file(GLOB_RECURSE  source_files src/* ccwrap/* pywrap/*)

# build test codes
enable_testing()
add_subdirectory(examples)

# define target library
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_library(apfel SHARED ${source_files})

# define libraries to be linked
target_link_libraries(apfel ${LHAPDF_LIBRARIES})

# install
install(TARGETS apfel
           RUNTIME DESTINATION bin
           LIBRARY DESTINATION lib)
install(DIRECTORY include/APFEL DESTINATION include)
install(FILES ${PROJECT_SOURCE_DIR}/bin/apfel-config
${CMAKE_CURRENT_BINARY_DIR}/examples/CheckAPFEL ${CMAKE_CURRENT_BINARY_DIR}/examples/ListFunctions DESTINATION bin
PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
