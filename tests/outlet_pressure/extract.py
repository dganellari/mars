"""Compile the actual boundary kernels on the host, without emulating a GPU solve."""

import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text()
output = [
    "#pragma once\n#include <cstdint>\n#define __global__\n"
    "struct Index { int x=0; };\ninline Index blockIdx, threadIdx;\ninline Index blockDim{1};\n"
]
for name in (
    "enforceBcMatrixKernel",
    "enforceBcColMatrixKernelByNode",
    "enforcePressureBcRhsLiftKernel",
):
    start = source.index("__global__ void " + name + "(")
    start = source.rindex("template<", 0, start)
    end = source.index("{", start) + 1
    depth = 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    output.append(source[start:end] + "\n")
pathlib.Path(sys.argv[2]).write_text("\n".join(output))
