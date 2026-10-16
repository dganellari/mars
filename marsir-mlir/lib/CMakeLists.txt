add_mlir_dialect_library(MLIRMir
  MirDialect.cpp
  MirOps.cpp
  LowerMirToLinalg.cpp
  HoistTransferPairs.cpp
  LowerCopies.cpp
  WarpDistribute.cpp
  WarpWrap.cpp
  RegisterResidentChain.cpp
  BatchElements.cpp
  GpuWrap.cpp
  ForwardTransfers.cpp
  EmulateWarp.cpp
  WorkgroupBuffers.cpp
  DistributeFills.cpp
  WarpBarriers.cpp
  HoistInvariantReads.cpp
  UnrollLoops.cpp
  ForwardOwned.cpp
  LowerCellwise.cpp
  WarpPerElement.cpp
  NoaliasDescriptors.cpp
  CellwiseSchedule.cpp
  FusePhases.cpp

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
  MLIRLinalgTransforms
  MLIRTensorDialect
  MLIRArithDialect
  MLIRSCFDialect
  MLIRSCFUtils
  MLIRVectorDialect
  MLIRMemRefDialect
  MLIRGPUDialect
  MLIRLLVMDialect
  MLIRNVGPUDialect
  MLIRFuncDialect
  )

# LowerCellwise.cpp takes the reference-hex numbering from MARS itself
# (backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp).
foreach(t MLIRMir obj.MLIRMir)
  if(TARGET ${t})
    target_include_directories(${t} PRIVATE ${PROJECT_SOURCE_DIR}/..)
  endif()
endforeach()
