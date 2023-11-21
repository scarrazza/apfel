![alt text](https://github.com/scarrazza/apfel/raw/master/resources/logoapfel.png "Logo APFEL")

# APFEL: A PDF Evolution Library

Visit: http://apfel.hepforge.org and http://apfel.mi.infn.it/

APFEL is a library able to perform DGLAP evolution up to NNLO in QCD
and to NLO in QED, both with pole and MSbar masses. The coupled DGLAP
QCD+QED evolution equations are solved in x-space by means of higher
order interpolations and Runge-Kutta techniques.

The APFEL library is accessible also through the APFEL Web
interface. APFEL Web provides an online web-application which
integrates several HEP softwares providing a complete suite plotting
tools for PDFs and many related quantities (http://apfel.mi.infn.it/).

## Download

You can obtain APFEL directly from the GitHub repository:

https://github.com/scarrazza/apfel/releases

For the last development version you can clone the master code:

```Shell
git clone https://github.com/scarrazza/apfel.git
```

For the latest tag:

```Shell
git tag -l
git checkout tags/tag_name
```

## Installation

### Dependencies

- Fortran and C++ compillers.
- (Optional) `CMake` (https://cmake.org/) > 3.16 for installation with `CMake` and any build system -- `make`, `ninja`, etc.
- (Optional) `LHAPDF` library, for support PDFs in of Les Houches format.
- (Optional) `CPython` headers and `SWIG`, for compilation of Python bindings.
- (Optional) Internet connection if the installation of PDFs was requested.

### From source

Checkout the code and compile the code by using the
following procedure:

```Shell
cd apfel
cmake -S . -B BUILD  <extra flags>
cmake --build BUILD
cmake --install BUILD
```

and optionally, if the testing was enabled

```Shell
ctest --test-dir BUILD
```

Remember to export the location of the libraries into the `LD_LIBRARY_PATH` or `DYLD_LIBRARY_PATH`.

The extra flags might be:

- generic CMake flags, e.g. `-DCMAKE_INSTALL_PREFIX=/my/home/dir`, `-DCMAKE_Fortran_COMPILER=ifort`, `-DCMAKE_Fortran_FLAGS="-O2 -g"`, etc.
- flags pointing to the dependencies, `-DPython_DIR=/where/the/derised/python/is`
- flags that steer the compilation. In the current version there are the following flags:
  - `-DAPFEL_ENABLE_PYTHON=ON|OFF`, default=`ON` Enables building of python bindings. Requires `SWIG` and Python headers if enabled. Tested only with CPython.
  - `-DAPFEL_ENABLE_TESTS=ON|OFF`, default=`ON` Enables testing
  - `-DAPFEL_ENABLE_LHAPDF=ON|OFF`, default=`ON` Enables compilation with `LHAPDF`. Requires `LHAPDF` if enabled.
  - `-DAPFEL_DOWNLOAD_PDFS=ON|OFF`, default=`ON` Download LHAPDF sets for tests. Makes sense only when the testing and compilation with LHAPDF are enabled. Requires internet connection if enabled.
  - `-DAPFEL_Python_SITEARCH=/path/to/install/python/modules|autoprefix`, default is the python system location. If "autoprefix" is used, the modules will be installed inside
    `CMAKE_INSTALL_PREFIX`. Makes sense only when the building of python bindings is enabled. Do not forget to add the location of modules to the `PYTHONPATH`.
    The `CMake` installation also provides the `CMake` config files for APFEL, therefore it if possible to
    do in the `CMakeLists.txt` of the dependant projects:

```Shell
find_package(apfel)
....
target_link_libraries(mytarget PRIVATE APFEL::APFEL APFEL::APFELevol)

```

#### Known issues

It is recommended to avoid using the source directory for the builds, i.e.

```Shell
cmake -S . -B . # do not do this
cmake -S . # do not do this
```

## References

- V. Bertone, S. Carrazza, J. Rojo, _APFEL: A PDF Evolution Library with QED corrections_, [arXiv:1310.1394](http://arxiv.org/abs/arXiv:1310.1394).
- S. Carrazza, A. Ferrara, D. Palazzo, J. Rojo, _APFEL Web: a web-based application for the graphical visualization of parton distribution functions_, [arXiv:1410.5456](http://arxiv.org/abs/1410.5456).

## Contact Information

Maintainers: Valerio Bertone, Stefano Carrazza

Homepage: http://apfel.hepforge.org/
