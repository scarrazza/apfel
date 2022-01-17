# BASIC DEFINITIONS ========================================================================

# define minimum version of cmake
cmake_minimum_required (VERSION 3.5)

cmake_policy(SET CMP0086 OLD)
cmake_policy(SET CMP0078 OLD)

# define project name, version and its languages
project(apfel CXX Fortran)
set(apfel_VERSION 3.0.7)

# define c++ standard and issue all the warning demanded by this standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -cpp")
set(CMAKE_SHARED_LINKER_FLAGS -w)

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# CONFIG SCRIPT =========================================================================

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
set(PACKAGE_VERSION "${apfel_VERSION}")
configure_file("${PROJECT_SOURCE_DIR}/bin/apfel-config.in" "${PROJECT_SOURCE_DIR}/bin/apfel-config")

# Fortran wrappers =========================================================================

configure_file("${PROJECT_SOURCE_DIR}/include/APFEL/FortranWrappers.h.cmake" "${PROJECT_SOURCE_DIR}/include/APFEL/FortranWrappers.h")

# FINALIZE ==================================================================================

# LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config)
if (LHAPDF_CONFIG)
  exec_program(${LHAPDF_CONFIG}
    ARGS --cxxflags
    OUTPUT_VARIABLE LHAPDF_CXX_FLAGS
  )
  set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${LHAPDF_CONFIG}
    ARGS --ldflags
    OUTPUT_VARIABLE LHAPDF_LIBRARIES
  )
  set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARIES} CACHE STRING INTERNAL)
else(LHAPDF_CONFIG)
  add_compile_definitions(NOLHAPDF)
  message("LHAPDF not found!")
endif(LHAPDF_CONFIG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS} -Wall -fPIC -Wunused")

# generate list of source files
file(GLOB_RECURSE  source_files src/* ccwrap/*)

# Add swig subdirectory
add_subdirectory(pywrap)

add_library(apfel SHARED ${source_files})

# include directory
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

# define libraries to be linked
target_link_libraries(apfel ${LHAPDF_LIBRARIES})

# build test codes
enable_testing()
add_subdirectory(examples)

# define target library
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# install
install(TARGETS apfel
           RUNTIME DESTINATION bin
           LIBRARY DESTINATION lib)
install(DIRECTORY include/APFEL DESTINATION include)
install(FILES ${PROJECT_SOURCE_DIR}/bin/apfel-config
${CMAKE_CURRENT_BINARY_DIR}/examples/CheckAPFEL ${CMAKE_CURRENT_BINARY_DIR}/examples/ListFunctions DESTINATION bin
PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

execute_process(COMMAND python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
install(TARGETS _apfel DESTINATION ${PYTHON_SITE_PACKAGES})
install(FILES ${CMAKE_BINARY_DIR}/pywrap/apfel.py DESTINATION ${PYTHON_SITE_PACKAGES})
