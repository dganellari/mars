#pragma once

// cuSPARSE and cuBLAS for CUDA and HIP builds. On AMD the generic SpMV and the BLAS calls the solvers use
// map to hipSPARSE and hipBLAS, which mirror them, so the solver code calls the CUDA names on both.
#include "mars_gpu_runtime.hpp"

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)

#include <hipblas/hipblas.h>
#include <hipsparse/hipsparse.h>

#define CUDA_R_32F HIP_R_32F
#define CUDA_R_64F HIP_R_64F

#define CUSPARSE_INDEX_32I HIPSPARSE_INDEX_32I
#define CUSPARSE_INDEX_BASE_ZERO HIPSPARSE_INDEX_BASE_ZERO
#define CUSPARSE_OPERATION_NON_TRANSPOSE HIPSPARSE_OPERATION_NON_TRANSPOSE
#define CUSPARSE_SPMV_ALG_DEFAULT HIPSPARSE_SPMV_ALG_DEFAULT
#define cusparseCreate hipsparseCreate
#define cusparseCreateCsr hipsparseCreateCsr
#define cusparseCreateDnVec hipsparseCreateDnVec
#define cusparseDestroy hipsparseDestroy
#define cusparseDestroyDnVec hipsparseDestroyDnVec
#define cusparseDnVecSetValues hipsparseDnVecSetValues
#define cusparseDestroySpMat hipsparseDestroySpMat
#define cusparseDnVecDescr_t hipsparseDnVecDescr_t
#define cusparseHandle_t hipsparseHandle_t
#define cusparseSpMatDescr_t hipsparseSpMatDescr_t
#define cusparseSpMV hipsparseSpMV
#define cusparseSpMV_bufferSize hipsparseSpMV_bufferSize

#define CUBLAS_OP_N HIPBLAS_OP_N
#define CUBLAS_OP_T HIPBLAS_OP_T
#define CUBLAS_POINTER_MODE_DEVICE HIPBLAS_POINTER_MODE_DEVICE
#define CUBLAS_POINTER_MODE_HOST HIPBLAS_POINTER_MODE_HOST
#define cublasCreate hipblasCreate
#define cublasDaxpy hipblasDaxpy
#define cublasDcopy hipblasDcopy
#define cublasDdot hipblasDdot
#define cublasDgemm hipblasDgemm
#define cublasDestroy hipblasDestroy
#define cublasDgemv hipblasDgemv
#define cublasDscal hipblasDscal
#define cublasHandle_t hipblasHandle_t
#define cublasSaxpy hipblasSaxpy
#define cublasScopy hipblasScopy
#define cublasSdot hipblasSdot
#define cublasSetPointerMode hipblasSetPointerMode
#define cublasSetStream hipblasSetStream
#define cublasSgemv hipblasSgemv
#define cublasSscal hipblasSscal

#else

#include <cublas_v2.h>
#include <cusparse.h>

#endif
