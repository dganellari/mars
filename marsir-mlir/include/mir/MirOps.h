#ifndef MIR_MIROPS_H
#define MIR_MIROPS_H

#include "mlir/Bytecode/BytecodeOpInterface.h"
#include "mlir/IR/BuiltinTypes.h"
#include "mlir/IR/Dialect.h"
#include "mlir/IR/OpDefinition.h"
#include "mlir/Interfaces/SideEffectInterfaces.h"

// Generated type declarations (the !mir.block type).
#define GET_TYPEDEF_CLASSES
#include "mir/MirOpsTypes.h.inc"

// Generated op declarations.
#define GET_OP_CLASSES
#include "mir/MirOps.h.inc"

#endif // MIR_MIROPS_H
