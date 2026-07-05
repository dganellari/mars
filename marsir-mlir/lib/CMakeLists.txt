add_mlir_dialect_library(MLIRMir
  MirDialect.cpp
  MirOps.cpp
  LowerMirToLinalg.cpp
  HoistTransferPairs.cpp
  LowerCopies.cpp
  WarpDistribute.cpp
  WarpWrap.cpp

  ADDITIONAL_HEADER_DIRS
  ${PROJECT_SOURCE_DIR}/include/mir

  DEPENDS
  MLIRMirOpsIncGen

  LINK_LIBS PUBLIC
  MLIRIR
  MLIRSupport
  MLIRPass
  MLIRTransforms
  MLIRTransformUtils
  MLIRLinalgDialect
  MLIRTensorDialect
  MLIRArithDialect
  MLIRSCFDialect
  MLIRVectorDialect
  MLIRMemRefDialect
  MLIRGPUDialect
  MLIRFuncDialect
  )
