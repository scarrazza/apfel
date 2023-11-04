# Set up swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

message(STATUS "APFEL python: SWIG version ${SWIG_VERSION} found in ${SWIG_EXECUTABLE}.")


# Include python
find_package(Python  COMPONENTS Interpreter Development REQUIRED)
message(STATUS "APFEL python: Python version ${Python_VERSION_MAJOR}.${Python_VERSION_MINOR} found in ${Python_EXECUTABLE}.")

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON  OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include/)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../include/)

# Add swig module
# Please note that the line below will generate apfelPYTHON_wrap.cxx file in the CMAKE_CURRENT_BINARY_DIR
swig_add_library(apfel LANGUAGE python SOURCES ${PROJECT_NAME}.i)
swig_link_libraries(apfel APFEL Python::Module)

# Write meta-info file
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/apfel-${apfel_VERSION}-py${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}.egg-info
"Metadata-Version: 2.1\n\
Name: APFEL\n\
Version: ${apfel_VERSION}\n\
Summary: A PDF Evolution Library\n\
Home-page: https://github.com/scarazza/apfel\n\
Author: ['Valerio Bertone']\n\
Author-email: valerio.bertone@cern.ch\n\
License: GPL\n\
Keywords: PDFs,FFs,DGLAP,DIS,structure functions\n\
\n\
This is a simple SWIG wrapper on the main steering interface of\n\
the APFEL library. It is used to create, query and use the APFEL data from\n\
a Python program.")

if (NOT APFEL_Python_SITEARCH)
  set(APFEL_Python_SITEARCH ${Python_SITEARCH})
  message(STATUS "APFEL python: WARNING: The installation path of the python modules is APFEL_Python_SITEARCH=${APFEL_Python_SITEARCH}.")
  message(STATUS "APFEL python: WARNING: The installation path of the python modules is outside of the global instalation path CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}.")
  message(STATUS "APFEL python: WARNING: The installation path of the python modules might be not writable by the current user.")
  message(STATUS "APFEL python: WARNING: You can use the APFEL_Python_SITEARCH variable to set the desired installation path for the Python modules.")
  message(STATUS "APFEL python: WARNING: Alternatively you can use -DAPFEL_Python_SITEARCH=autoprefix to put the python modules in the CMAKE_INSTALL_PREFIX directory.")
elseif ( "${APFEL_Python_SITEARCH}" STREQUAL "autoprefix")
  set(APFEL_Python_SITEARCH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/")
  message(STATUS "APFEL python: WARNING: The installation path of the python modules is APFEL_Python_SITEARCH=${APFEL_Python_SITEARCH} is inside the CMAKE_INSTALL_PREFIX.")
  message(STATUS "APFEL python: WARNING: Don't forget to add this location to your PYTHONPATH.")
else()
  message(STATUS "APFEL python: The installation path of the python modules is APFEL_Python_SITEARCH=${APFEL_Python_SITEARCH} is set by the user.")
  message(STATUS "APFEL python: Don't forget to add this location to your PYTHONPATH.")
  if (NOT IS_ABSOLUTE ${APFEL_Python_SITEARCH})
    message(STATUS "APFEL python: WARNING: The path APFEL_Python_SITEARCH=${APFEL_Python_SITEARCH} is not absolute! It is recommented to use absolute path!")
  endif()
endif()
# Do manual instalation and don't rely on the python tools
# This renders the setup.py.in obsolete.
install(TARGETS _apfel DESTINATION ${APFEL_Python_SITEARCH}/apfel)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/apfel-${apfel_VERSION}-py${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}.egg-info DESTINATION ${APFEL_Python_SITEARCH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/apfel.py DESTINATION ${APFEL_Python_SITEARCH}/apfel)
