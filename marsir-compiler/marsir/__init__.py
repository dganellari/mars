"""marsir -- MARS FEM operator IR / codegen (Stage 1: source emitter, no LLVM).

Pipeline: frontend.parse_spec -> ir.synthesize -> backends.{cuda,host_cpp}.emit
"""
from .frontend import parse_spec, parse_spec_file
from .ir import synthesize, dump, Operator, ElementApply, Block, block_for

__all__ = ["parse_spec", "parse_spec_file", "synthesize", "dump",
           "Operator", "ElementApply", "Block", "block_for"]
