# Set up swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# Include python
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
include_directories(../include/)

# Add swig module
swig_add_library(apfel LANGUAGE python SOURCES ${PROJECT_NAME}.i)
swig_link_libraries(apfel ${PROJECT_NAME} ${PYTHON_LIBRARIES})

# Files to install with Python
set(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/apfel.py
        ${CMAKE_CURRENT_BINARY_DIR}/_apfel.so)

# Configure setup.py and copy to output directory
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})

# Install target to call setup.py
add_custom_target(install-python
        DEPENDS _apfel
        COMMAND python ${SETUP_PY_OUT} install)
