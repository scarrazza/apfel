# BASIC DEFINITIONS ========================================================================

# define minimum version of cmake
cmake_minimum_required (VERSION 3.15)

cmake_policy(SET CMP0086 OLD)
cmake_policy(SET CMP0078 OLD)

# define project name, version and its languages
project(apfel LANGUAGES CXX Fortran)
set(apfel_VERSION 3.0.7)
set(CMAKE_VERBOSE_MAKEFILE OFF)
#This module respects HFS, e.g. defines lib or lib64 when it is needed.
include("GNUInstallDirs")

option(APFEL_ENABLE_PYTHON      "Enables building of python bindings" ON)
option(APFEL_ENABLE_TESTS       "Enables testing" ON)
option(APFEL_ENABLE_LHAPDF      "Enables LHAPDF" ON)
option(APFEL_DOWNLOAD_PDFS      "Download LHAPDF sets" ON)

set(CMAKE_BUILD_TYPE Release)
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wall>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wextra>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-unused-parameter>")

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# CONFIG SCRIPT =========================================================================

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
set(libdir "${prefix}/${CMAKE_INSTALL_LIBDIR}")
set(PACKAGE_VERSION "${apfel_VERSION}")
configure_file("${PROJECT_SOURCE_DIR}/bin/apfel-config.in" "${PROJECT_BINARY_DIR}/bin/apfel-config")

# Fortran wrappers =========================================================================

configure_file("${PROJECT_SOURCE_DIR}/include/APFEL/FortranWrappers.h.cmake" "${PROJECT_BINARY_DIR}/include/APFEL/FortranWrappers.h")

# FINALIZE ==================================================================================

# LHAPDF
if (APFEL_ENABLE_LHAPDF)
find_program(LHAPDF_CONFIG lhapdf-config)
endif()
if (LHAPDF_CONFIG AND APFEL_ENABLE_LHAPDF)
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
else()
  add_compile_definitions(NOLHAPDF)
  message("LHAPDF is disabled of not found!")
endif()

# generate list of source files
file(GLOB_RECURSE  source_files_cxx src/*/*.cc ccwrap/*.cc)
file(GLOB_RECURSE  source_files_fortran src/*/*.f ccwrap/*.f)
set_source_files_properties( ${source_files_fortran}  PROPERTIES Fortran_PREPROCESS ON )

# Add swig subdirectory
if (APFEL_ENABLE_PYTHON)
  add_subdirectory(pywrap)
endif()

add_library(APFEL SHARED ${source_files_cxx} ${source_files_fortran})

# include directory
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

# define libraries to be linked
target_link_libraries(APFEL ${LHAPDF_LIBRARIES})
target_compile_features(APFEL PRIVATE cxx_std_11)
set_target_properties(APFEL PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    )

# build test codes
if (APFEL_ENABLE_TESTS)
 enable_testing()
 add_subdirectory(examples)
endif()

# install
install(TARGETS APFEL DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY include/APFEL DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")
install(FILES ${PROJECT_BINARY_DIR}/include/APFEL/FortranWrappers.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/APFEL/ )
install(FILES ${PROJECT_BINARY_DIR}/bin/apfel-config
${CMAKE_CURRENT_BINARY_DIR}/examples/CheckAPFEL ${CMAKE_CURRENT_BINARY_DIR}/examples/ListFunctions DESTINATION ${CMAKE_INSTALL_BINDIR}
PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

