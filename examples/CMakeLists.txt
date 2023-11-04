# build tests

file(GLOB testcodes *.cc *.f)
foreach(testsource ${testcodes})
  GET_FILENAME_COMPONENT(filename ${testsource} NAME_WE)
  add_executable(${filename} ${testsource})
  target_link_libraries(${filename} APFEL)
  if (${filename} STREQUAL "Timing" )
    add_test(NAME ${filename}  COMMAND sh -c "${CMAKE_CURRENT_BINARY_DIR}/${filename} 2 QCD" )
  elseif (${filename} STREQUAL "TabulationStep" )
    #add_test(NAME ${filename}  COMMAND sh -c "${CMAKE_CURRENT_BINARY_DIR}/${filename}  <${CMAKE_CURRENT_SOURCE_DIR}/input.txt" )
  elseif (${filename} STREQUAL "TabulationStepCxx" )
    #add_test(NAME ${filename}  COMMAND sh -c "${CMAKE_CURRENT_BINARY_DIR}/${filename}  <${CMAKE_CURRENT_SOURCE_DIR}/input.txt" )
  else()
    add_test(NAME ${filename}  COMMAND sh -c "${CMAKE_CURRENT_BINARY_DIR}/${filename} <${CMAKE_CURRENT_SOURCE_DIR}/input.txt" )
  endif()
endforeach()


if (APFEL_DOWNLOAD_PDFS AND APFEL_ENABLE_LHAPDF)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PDFS)
  file(DOWNLOAD https://lhapdfsets.web.cern.ch/current/NNPDF23_nlo_as_0118.tar.gz ${CMAKE_CURRENT_BINARY_DIR}/PDFS/NNPDF23_nlo_as_0118.tar.gz)
  file(DOWNLOAD https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_as_0118.tar.gz ${CMAKE_CURRENT_BINARY_DIR}/PDFS/NNPDF31_nnlo_as_0118.tar.gz)
  add_custom_target( UNZIPPDFS ALL)
  add_custom_command(TARGET UNZIPPDFS PRE_BUILD 
         COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/PDFS/NNPDF23_nlo_as_0118
         COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/PDFS/NNPDF31_nnlo_as_0118
         COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_BINARY_DIR}/PDFS/NNPDF23_nlo_as_0118.tar.gz
         COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_BINARY_DIR}/PDFS/NNPDF31_nnlo_as_0118.tar.gz
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PDFS
         DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/PDFS/NNPDF23_nlo_as_0118.tar.gz ${CMAKE_CURRENT_BINARY_DIR}/PDFS/NNPDF31_nnlo_as_0118.tar.gz
         COMMENT "Unpacking PDFs"
         VERBATIM)
  foreach( t LHgridDerivativeProduction LHgridProduction LHgridDerivativeProductionCxx)
    set_tests_properties( ${t} PROPERTIES ENVIRONMENT "LHAPDF_DATA_PATH=${CMAKE_CURRENT_BINARY_DIR}/PDFS" DEPENDS UNZIPPDFS)
  endforeach()

endif()
