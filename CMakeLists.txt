# BASIC DEFINITIONS ====================================================

# define minimum version of cmake
cmake_minimum_required (VERSION 3.15)

cmake_policy(SET CMP0086 OLD)
cmake_policy(SET CMP0078 OLD)

# define project name, version and its languages
project(apfel LANGUAGES C CXX Fortran)
set(apfel_VERSION 3.1.1)
set(PROJECT_VERSION 3.1.1)
SET(APFEL_VERSION_MAJOR  3)
SET(APFEL_VERSION_MINOR  1)
SET(APFEL_VERSION_PATCH  1)
add_compile_definitions(APFEL_VERSION=${apfel_VERSION})
set(CMAKE_VERBOSE_MAKEFILE OFF)
include("GNUInstallDirs")

option(APFEL_ENABLE_PYTHON      "Enables building of python bindings" ON)
option(APFEL_ENABLE_TESTS       "Enables testing" ON)
option(APFEL_ENABLE_LHAPDF      "Enables LHAPDF" ON)
option(APFEL_DOWNLOAD_PDFS      "Download LHAPDF sets" ON)

set(CMAKE_BUILD_TYPE Release)
# Flags for GNU compilers ==============================================
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wall>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wextra>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-unused-parameter>")

# RPATH  +++++++++++++++++==============================================
if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()


# Output of config info   ==============================================
message(STATUS "APFEL: CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}")
message(STATUS "APFEL: CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
message(STATUS "APFEL: CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "APFEL: APFEL_ENABLE_PYTHON=${APFEL_ENABLE_PYTHON}")
message(STATUS "APFEL: APFEL_ENABLE_TESTS=${APFEL_ENABLE_TESTS}")
message(STATUS "APFEL: APFEL_ENABLE_LHAPDF=${APFEL_ENABLE_LHAPDF}")
message(STATUS "APFEL: APFEL_DOWNLOAD_PDFS=${APFEL_DOWNLOAD_PDFS}")

# CONFIG SCRIPT ========================================================
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(includedir "${CMAKE_INSTALL_FULL_INCLUDEDIR}")
set(libdir "${CMAKE_INSTALL_FULL_LIBDIR}")
set(PACKAGE_VERSION "${apfel_VERSION}")
configure_file("${PROJECT_SOURCE_DIR}/bin/apfel-config.in" "${PROJECT_BINARY_DIR}/bin/apfel-config")
configure_file("${PROJECT_SOURCE_DIR}/bin/apfel.in" "${PROJECT_BINARY_DIR}/bin/apfel")

# Fortran wrappers =====================================================

include(FortranCInterface)
message(STATUS "APFEL: FortranCInterface FortranCInterface_GLOBAL_CASE=${FortranCInterface_GLOBAL_CASE},  FortranCInterface_GLOBAL__CASE=${FortranCInterface_GLOBAL__CASE}")
if ("${FortranCInterface_GLOBAL_CASE}" STREQUAL "UPPER")
 set(N "NAME")
else()
 set(N "name")
endif()
if ("${FortranCInterface_GLOBAL__CASE}" STREQUAL "UPPER")
 set(N_ "NAME")
else()
 set(N_ "name")
endif()
FILE(WRITE   "${PROJECT_BINARY_DIR}/include/APFEL/FortranWrappers.h"
"#ifndef APFEL_FORTRANWRAPPERS_H\n\
#define APFEL_FORTRANWRAPPERS_H\n\
\n\
/* Mangling for Fortran global symbols without underscores. */\n\
#define FC_FUNC(name,NAME) ${FortranCInterface_GLOBAL_PREFIX}${N}##${FortranCInterface_GLOBAL_SUFFIX} \n\
\n\
/* Mangling for Fortran global symbols with underscores. */\n\
#define FC_FUNC_(name,NAME) ${FortranCInterface_GLOBAL__PREFIX}${N_}##${FortranCInterface_GLOBAL__SUFFIX} \n\
#endif\n")

# LHAPDF ===============================================================
if (APFEL_ENABLE_LHAPDF)
  find_program(LHAPDF_CONFIG lhapdf-config)
  if (LHAPDF_CONFIG)
    message(STATUS "APFEL: lhapdf-config found in  ${LHAPDF_CONFIG}.")
  endif()
endif()
if (LHAPDF_CONFIG AND APFEL_ENABLE_LHAPDF)
  exec_program(${LHAPDF_CONFIG}
    ARGS --cxxflags
    OUTPUT_VARIABLE LHAPDF_CXX_FLAGS
  )
  set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${LHAPDF_CONFIG}
    ARGS --libdir
    OUTPUT_VARIABLE LHAPDF_LIBDIR
  )
  find_library(LHAPDF_LIBRARY NAMES LHAPDF PATHS ${LHAPDF_LIBDIR})
  set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARY} CACHE STRING INTERNAL)

  exec_program(${LHAPDF_CONFIG}
    ARGS --incdir
    OUTPUT_VARIABLE LHAPDF_INCLUDE_DIRS
  )
  set(LHAPDF_INCLUDE_DIRS ${LHAPDF_INCLUDE_DIRS} CACHE STRING INTERNAL)
  message(STATUS "APFEL: LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES}")
else()
  add_compile_definitions(NOLHAPDF)
  message("LHAPDF is disabled of not found!")
endif()

# List of source files =================================================

set(Fortran_PREPROCESS ON )


set(INITCXX  Evolution/apfelversion.cc Evolution/lhapdfwrapper.cc)

set(INIT      Evolution/WelcomeMessage.f Evolution/initParameters.f Evolution/initGrid.f Evolution/initGridAlpha.f 
    Evolution/initIntegralsQCD.f Evolution/initIntegralsQCDRes.f Evolution/initIntegralsQED.f 
    Evolution/initIntegralsMatching.f Evolution/initPDFs.f Evolution/InitializeAPFEL.f Evolution/EvolveAPFEL.f 
    Evolution/DeriveAPFEL.f Evolution/initHELL.f Evolution/TruncatedEvolveAPFEL.f Evolution/ExponentiatedEvolveAPFEL.f 
    Evolution/ReportParameters.f Evolution/CheckAPFEL.f Evolution/ComputeHeavyQuarkThresholds.f 
    Evolution/CachePDFsAPFEL.f)

set(SETTINGS  Evolution/SetAlphaQCDRef.f Evolution/SetAlphaQEDRef.f Evolution/SetQLimits.f 
    Evolution/SetPerturbativeOrder.f Evolution/SetFFNS.f Evolution/SetGridParameters.f Evolution/SetExternalGrid.f 
    Evolution/SetMSbarMasses.f Evolution/SetNumberOfGrids.f Evolution/SetPDFSet.f Evolution/SetPoleMasses.f 
    Evolution/SetRenFacRatio.f Evolution/SetTheory.f Evolution/SetVFNS.f Evolution/SetMaxFlavourPDFs.f 
    Evolution/SetMaxFlavourAlpha.f Evolution/SetReplica.f Evolution/AlphaQCD.f Evolution/AlphaQED.f Evolution/HeavyQuarkMass.f 
    Evolution/LHAPDFgrid.f Evolution/CleanUp.f Evolution/EnableWelcomeMessage.f Evolution/EnableEvolutionOperator.f 
    Evolution/ExternalEvolutionOperator.f Evolution/LockGrids.f Evolution/SetTimeLikeEvolution.f Evolution/SetAlphaEvolution.f 
    Evolution/SetLambdaQCDRef.f Evolution/SetPDFEvolution.f Evolution/LHAPDFgridDerivative.f 
    Evolution/SetSmallxResummation.f Evolution/SetFastEvolution.f Evolution/EnableMassRunning.f Evolution/ListFunctions.f 
    Evolution/GetPerturbativeOrder.f Evolution/SetTauMass.f Evolution/EnableLeptonEvolution.f Evolution/SetLHgridParameters.f 
    Evolution/SetEpsilonTruncation.f Evolution/GetMuF.f Evolution/SetPolarizedEvolution.f Evolution/SetMassMatchingScales.f 
    Evolution/SetMassScaleReference.f Evolution/SetQGridParameters.f Evolution/GetThreshold.f Evolution/EnableNLOQEDCorrections.f 
    Evolution/ExternalSplittingFunctions.f Evolution/GetMaxFlavourPDFs.f Evolution/GetMaxFlavourAlpha.f Evolution/GetLambdaQCD.f)

set(INTER  Evolution/interpolants.f)

set(SPLITTING_FUNCTIONS  Evolution/SplittingFunctions.f)

set(COUPLINGS  Evolution/a_QCD.f Evolution/a_QED.f Evolution/MSbarmass.f Evolution/lambda.f)

set(EVOLUTION  Evolution/EvolutionQCD.f Evolution/EvolutionOperatorsQCD.f Evolution/IdentityOperatorsQCD.f 
    Evolution/EqualOperatorsQCDnf.f Evolution/integrandsQCD.f Evolution/integralsQCD.f Evolution/RSLintegralsQCD.f 
    Evolution/odeintnsQCD.f Evolution/odeintsgQCD.f Evolution/MatchingConditions.f Evolution/integrandsMatching.f 
    Evolution/integralsMatching.f Evolution/RSLintegralsMatching.f Evolution/MatchPDFs.f Evolution/JoinOperatorsQCD.f 
    Evolution/EvolutionUnified.f Evolution/EvolutionOperatorsUnified.f Evolution/IdentityOperatorsUnified.f 
    Evolution/EqualOperatorsUnifiednf.f Evolution/EvolveUni.f Evolution/odeintnsUnified.f Evolution/odeintsgUnified.f 
    Evolution/DerivativeOperatorsQCD.f Evolution/odeintnsQCDf.f Evolution/odeintsgQCDf.f Evolution/odeintnsUnifiedf.f 
    Evolution/odeintsgUnifiedf.f Evolution/integrandsQED.f Evolution/integralsQED.f Evolution/RSLintegralsQED.f 
    Evolution/JoinOperatorsUni.f Evolution/ExternalEvolutionMatrixEv2Ev.f Evolution/ExternalEvolutionMatrixEv2Ph.f 
    Evolution/ExternalEvolutionMatrixPh2Ph.f)

set(PDF  Evolution/PDFphys2evQCD.f Evolution/PDFevQCD2phys.f Evolution/EvolvePDFs.f Evolution/EvolveQCD.f Evolution/toyLHPDFs.f 
    Evolution/xPDF.f Evolution/xgamma.f Evolution/NPDF.f Evolution/Ngamma.f Evolution/LUMI.f Evolution/JoinGrids.f 
    Evolution/PDFevUni2evQCD.f Evolution/PDFevQCD2evUni.f Evolution/PDFevUni2phys.f Evolution/PDFphys2evUni.f 
    Evolution/dxPDF.f Evolution/dxgamma.f Evolution/DeriveQCD.f Evolution/DerivePDFs.f Evolution/xGrid.f 
    Evolution/nIntervals.f Evolution/xLepton.f Evolution/xPDFxQ.f)

set(UTILITIES  core/dgauss.f core/wgplg.f core/ddilog.f core/gamma.f core/beta.f core/hplog.f)

set(DIS    DIS/InitializeAPFEL_DIS.f DIS/initParametersDIS.f DIS/SetMassScheme.f DIS/initIntegralsDIS.f 
    DIS/ZeroMassCoefficientFunctions.f DIS/RSLintegralsDIS.f DIS/ComputeStructureFunctionsAPFEL.f DIS/F2total.f 
    DIS/F2light.f DIS/F2charm.f DIS/F2bottom.f DIS/F2top.f DIS/FLtotal.f DIS/FLlight.f DIS/FLcharm.f 
    DIS/FLbottom.f DIS/FLtop.f DIS/F3total.f DIS/F3light.f DIS/F3charm.f DIS/F3bottom.f DIS/F3top.f 
    DIS/integrandsDIS.f DIS/SetProcessDIS.f DIS/ComputeChargesDIS.f DIS/SetPolarizationDIS.f DIS/SetProjectileDIS.f 
    DIS/SetTargetDIS.f DIS/MassiveCoefficientFunctions.f DIS/MassiveZeroCoefficientFunctions.f DIS/ComputeDISOperators.f 
    DIS/ConvolutePDFsWithDISOperators.f DIS/JoinDISOperators.f DIS/ConvoluteEvolutionWithDISOperators.f 
    DIS/ExternalDISOperator.f DIS/SetCKM.f DIS/SetGFermi.f DIS/SetProtonMass.f DIS/SetSin2ThetaW.f 
    DIS/SetWMass.f DIS/SetZMass.f DIS/FKSimulator.f DIS/SetFKObservable.f DIS/GetFKObservable.f DIS/FKObservables.f 
    DIS/EnableTargetMassCorrections.f DIS/GetCKM.f DIS/GetGFermi.f DIS/GetProtonMass.f DIS/GetSin2ThetaW.f 
    DIS/GetWMass.f DIS/GetZMass.f DIS/RSLintegralsSIA.f DIS/integrandsSIA.f DIS/initIntegralsSIA.f DIS/SelectCharge.f 
    DIS/SetRenQRatio.f DIS/SetFacQRatio.f DIS/GetSIATotalCrossSection.f DIS/IncludeScaleVariation.f 
    DIS/ExcludeScaleVariation.f DIS/ReportParametersDIS.f DIS/EnableDampingFONLL.f DIS/SetPropagatorCorrection.f 
    DIS/SetEWCouplings.f DIS/EnableDynamicalScaleVariations.f DIS/EnableIntrinsicCharm.f DIS/IncludeIntrinsicCharm.f 
    DIS/integrandsIC.f DIS/CacheStructureFunctionsAPFEL.f DIS/StructureFunctionxQ.f DIS/SetDampingPowerFONLL.f 
    DIS/hqcoef.f DIS/initIntegralsDISRes.f DIS/IncludeNLOQEDCorrections.f DIS/EnableSFNLOQEDCorrections.f 
    DIS/LHAPDFgridStructureFunctions.f DIS/SetScaleVariationProcedure.f DIS/SetNCComponent.f 
    DIS/RSLintegralspDIS.f DIS/integrandspDIS.f DIS/initIntegralspDIS.f 
    DIS/g4total.f DIS/g4light.f DIS/g4charm.f DIS/g4bottom.f DIS/g4top.f 
    DIS/g1total.f DIS/g1light.f DIS/g1charm.f DIS/g1bottom.f DIS/g1top.f 
    DIS/gLtotal.f DIS/gLlight.f DIS/gLcharm.f DIS/gLbottom.f DIS/gLtop.f)

set(FTDY  FTDY/src/ComputeFKTables.f FTDY/src/ComputeHardCrossSectionsDY.f FTDY/src/gauleg.f FTDY/src/initDYcouplings.f FTDY/src/initxGridDY.f 
    FTDY/src/readcDY.f FTDY/src/ReadDataFile.f FTDY/src/sigmafk_dy.f FTDY/src/wcDYlo.f FTDY/src/wcDYqg.f FTDY/src/wcDYqq.f FTDY/src/wcoeffDY.f 
    FTDY/src/writecDY.f FTDY/src/writeFK.f FTDY/src/xlineigen.f FTDY/src/PDFphys2ev.f)

set(HELL  HELL/src/hell-x.cc HELL/src/hell_wrapper.cc HELL/src/expansionSFs.cc HELL/src/math/special_functions.cc)


set(libAPFELCore_la_SOURCES  ${INITCXX} ${INIT} ${SETTINGS} ${INTER} ${SPLITTING_FUNCTIONS} ${COUPLINGS} ${EVOLUTION} ${PDF} ${UTILITIES} ${DIS} ${HELL} ${FTDY})
list(TRANSFORM libAPFELCore_la_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/")
set(libAPFELevol_la_SOURCES  ${INIT} ${SETTINGS} ${INTER} ${SPLITTING_FUNCTIONS} ${COUPLINGS} ${EVOLUTION} ${PDF} ${UTILITIES} ${HELL})
list(TRANSFORM libAPFELevol_la_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/")
set(libAPFELCCWrapEvol_la_SOURCES  ccwrap/APFELevol.cc )
set(libAPFELCCWrapObs_la_SOURCES  ccwrap/APFELobs.cc ccwrap/APFELdev.cc)

## Include directories #################################################

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)



# Add swig subdirectory ################################################
if (APFEL_ENABLE_PYTHON)
  add_subdirectory(pywrap)
endif()


# Shared libraries  ####################################################
add_library(APFEL SHARED ${libAPFELCCWrapObs_la_SOURCES} ${libAPFELCCWrapEvol_la_SOURCES} ${libAPFELCore_la_SOURCES})
target_include_directories(APFEL PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>  ${LHAPDF_INCLUDE_DIRS})
target_link_libraries(APFEL ${LHAPDF_LIBRARIES})
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  target_link_options(APFEL PRIVATE "LINKER:-z,defs")
  target_link_options(APFEL PRIVATE "LINKER:-z,now")
endif()
target_compile_features(APFEL PRIVATE cxx_std_11)
set_target_properties(APFEL PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    SOVERSION 0.0.0
    EXPORT_NAME APFEL::APFEL
    )
add_library(APFEL::APFEL ALIAS APFEL)

add_library(APFELevol SHARED ${libAPFELCCWrapEvol_la_SOURCES} ${libAPFELevol_la_SOURCES})
target_include_directories(APFELevol PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>  ${LHAPDF_INCLUDE_DIRS})
target_link_libraries(APFELevol ${LHAPDF_LIBRARIES} APFEL)
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  target_link_options(APFELevol PRIVATE "LINKER:-z,defs")
  target_link_options(APFELevol PRIVATE "LINKER:-z,now")
endif()
target_compile_features(APFELevol PRIVATE cxx_std_11)
set_target_properties(APFELevol PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    SOVERSION 0.0.0
    EXPORT_NAME APFEL::APFELevol
    )
add_library(APFEL::APFELevol ALIAS APFELevol)

# Static libraries  ####################################################
add_library(APFEL_static STATIC ${libAPFELCCWrapObs_la_SOURCES} ${libAPFELCCWrapEvol_la_SOURCES} ${libAPFELCore_la_SOURCES})
target_include_directories(APFEL_static PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>  ${LHAPDF_INCLUDE_DIRS})
target_link_libraries(APFEL_static ${LHAPDF_LIBRARIES})
target_compile_features(APFEL_static PRIVATE cxx_std_11)
set_target_properties(APFEL_static PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    OUTPUT_NAME APFEL
    EXPORT_NAME APFEL::APFEL_static
    )

add_library(APFELevol_static STATIC ${libAPFELCCWrapEvol_la_SOURCES} ${libAPFELevol_la_SOURCES})
target_include_directories(APFELevol_static PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>  ${LHAPDF_INCLUDE_DIRS})
target_link_libraries(APFELevol_static ${LHAPDF_LIBRARIES})
target_compile_features(APFELevol_static PRIVATE cxx_std_11)
set_target_properties(APFELevol_static PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    OUTPUT_NAME APFELevol
    EXPORT_NAME APFEL::APFELevol_static
    )



# Build test codes   ###################################################
if (APFEL_ENABLE_TESTS)
 enable_testing()
 add_subdirectory(examples)
endif()


# Cmake config files   #################################################
include(CMakePackageConfigHelpers)

set(CONFIG_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/APFEL/cmake)

configure_package_config_file(cmake/Templates/APFELConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/share/APFEL/cmake/APFELConfig.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/share/APFEL/cmake/APFELConfig-version.cmake COMPATIBILITY SameMajorVersion VERSION ${PROJECT_VERSION})


install(FILES ${CMAKE_CURRENT_BINARY_DIR}/share/APFEL/cmake/APFELConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/share/APFEL/cmake/APFELConfig-version.cmake
        DESTINATION ${CONFIG_INSTALL_DIR} COMPONENT devel)


# Install  #############################################################
## Libraries
install(TARGETS APFEL APFELevol APFEL_static APFELevol_static   EXPORT APFELTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
## Headers
install(DIRECTORY include/APFEL DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")
install(FILES ${PROJECT_BINARY_DIR}/include/APFEL/FortranWrappers.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/APFEL/ )
## Config scripts
install(FILES ${PROJECT_BINARY_DIR}/bin/apfel-config ${PROJECT_BINARY_DIR}/bin/apfel 
                                DESTINATION ${CMAKE_INSTALL_BINDIR} 
                                PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
## Readme
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/README.md DESTINATION  ${CMAKE_INSTALL_DOCDIR}/ COMPONENT doc)
## CMake config
install(EXPORT APFELTargets DESTINATION ${CONFIG_INSTALL_DIR} COMPONENT devel)


