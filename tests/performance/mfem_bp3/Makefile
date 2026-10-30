# Build against an MFEM install: make -f <this file> MFEM_DIR=<prefix holding share/mfem/config.mk>.
# MFEM's config.mk gives the compiler (nvcc for a CUDA build), flags and libraries. The
# binary goes to the current directory.
SRC_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include $(MFEM_DIR)/share/mfem/config.mk

mfem_bp3: $(SRC_DIR)mfem_bp3.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

clean:
	rm -f mfem_bp3

.PHONY: clean
