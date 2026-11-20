get_property(dialect_libs GLOBAL PROPERTY MLIR_DIALECT_LIBS)
get_property(conversion_libs GLOBAL PROPERTY MLIR_CONVERSION_LIBS)
get_property(extension_libs GLOBAL PROPERTY MLIR_EXTENSION_LIBS)

add_llvm_executable(mir-opt mir-opt.cpp)
llvm_update_compile_flags(mir-opt)

target_link_libraries(mir-opt PRIVATE
  ${dialect_libs}
  ${conversion_libs}
  ${extension_libs}
  MLIRMir
  MLIROptLib
  MLIRToLLVMIRTranslationRegistration
  MLIRIR
  MLIRSupport
  )

mlir_check_all_link_libraries(mir-opt)
