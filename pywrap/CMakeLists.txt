# Set up swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# Include python
find_package(Python  COMPONENTS Interpreter Development REQUIRED)

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON  OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include/)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../include/)

# Add swig module
swig_add_library(apfel LANGUAGE python SOURCES ${PROJECT_NAME}.i)
swig_link_libraries(apfel APFEL Python::Module)

# Files to install with Python
set(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/apfel.py
        ${CMAKE_CURRENT_BINARY_DIR}/_apfel.so)

# Configure setup.py and copy to output directory

set(srcdir ${CMAKE_CURRENT_BINARY_DIR})
set(top_srcdir ${PROJECT_SOURCE_DIR})
set(top_builddir ${PROJECT_BINARY_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in ${CMAKE_CURRENT_BINARY_DIR}/setup.py)


if (NOT APFEL_Python_SITEARCH)
 set(APFEL_Python_SITEARCH ${Python_SITEARCH})
endif()
if ( "${APFEL_Python_SITEARCH}" STREQUAL "prefix")
 set(APFEL_Python_SITEARCH "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/")
endif()
install(TARGETS _apfel DESTINATION ${APFEL_Python_SITEARCH})
install(FILES ${CMAKE_BINARY_DIR}/pywrap/apfel.py DESTINATION ${APFEL_Python_SITEARCH})


# Install target to call setup.py
add_custom_target(install-python
        DEPENDS _apfel
        COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py  install --prefix=${APFEL_Python_SITEARCH})