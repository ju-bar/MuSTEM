# Makefile for muSTEM
#
# === Usage ===	
# By Default, this Makefile builds the muSTEM application for CPU with single precision.
# You can customize the build by setting the `MODE` and `PRECISION` variables.
# The `MODE` variable can be set to `CPU` or `GPU`, and the `PRECISION` variable can be set to `single` or `double`.
#
# Take care to initialize the Compiler and CUDA environment variables correctly before running the Makefile.
#   source /usr/local/intel/oneapi/setvars.sh
#    or
#   source /opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin/nvhpcvars.sh
#  before running the Makefile (you can try first and run those if  it fails)
#
# To build for CPU:
#   make
# To build for GPU:
#   make MODE=GPU
# To build with double precision:
#   make MODE=CPU PRECISION=double
#  or
#   make MODE=GPU PRECISION=double
# To clean build artifacts:
#   make clean
# To rebuild everything:
#   make rebuild
#
# === muSTEM Makefile ===
# === Default Configuration ===
MODE      ?= CPU
PRECISION ?= single

# === Compiler & Toolchain Configuration ===
NVCOMP_BIN := /opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin
CUDA_ROOT  := /opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda

# === GPU Mode Configuration ===
ifeq ($(MODE),GPU)
  FC := $(NVCOMP_BIN)/nvfortran
  #FCFLAGS := -Mpreprocess -DGPU -Mbackslash -Mfree -mp -cuda -fast -O3
  FCFLAGS := -Mpreprocess -DGPU -Mbackslash -Mfree -mp -cuda -O2 -Kieee
  INCLUDES := -I$(CUDA_ROOT)/include
  LIBDIRS := -L$(CUDA_ROOT)/lib64
  LIBS := -lcufft -lrt -lpthread -ldl
  OBJDIR := x64/objgpu
  OUTNAME := muSTEM_GPU

else
# === CPU Mode Configuration ===
  ifndef MKLROOT
    $(error Please run 'source /usr/local/intel/oneapi/setvars.sh' before calling make)
  endif
  FC := ifx
  FCFLAGS := -fpp -qopenmp -O2 -heap-arrays -extend-source 132 -assume byterecl
  INCLUDES := -I$(MKLROOT)/include/fftw
  LIBDIRS := -L$(MKLROOT)/lib/intel64
  LIBS := -qmkl -liomp5
  OBJDIR := x64/objcpu
  OUTNAME := muSTEM
endif

# === Precision Flag ===
ifeq ($(PRECISION),double)
  FCFLAGS += -Ddouble_precision
else
  FCFLAGS += -Dsingle_precision
endif

# === Output Directory ===
OUTDIR := x64/Release

# === Source Files ===
GPU_MOD_SRC := \
	m_precision.f90 \
	mod_global_variables.f90 \
	m_string.f90 \
	m_user_input.f90 \
	GPU_routines/mod_cufft.f90 \
	GPU_routines/mod_cuda_setup.f90 \
	GPU_routines/mod_cuda_array_library.f90 \
	mod_CUFFT_wrapper.f90 \
	mod_output.f90 \
	m_numerical_tools.f90 \
	m_crystallography.f90 \
	m_plasmon.f90 \
	m_multislice.f90 \
	m_lens.f90 \
	m_tilt.f90 \
	quadpack.f90 \
	m_electron.f90 \
	m_absorption.f90 \
	GPU_routines/mod_cuda_potential.f90 \
	m_potential.f90 \
	mod_Hn0.f90 \
	GPU_routines/mod_cuda_ms.f90 \
	MS_utilities.f90

CPU_MOD_SRC := \
	m_precision.f90 \
	mod_global_variables.f90 \
	m_string.f90 \
	m_user_input.f90 \
	mod_CUFFT_wrapper.f90 \
	mod_output.f90 \
	m_numerical_tools.f90 \
	m_crystallography.f90 \
	m_plasmon.f90 \
	m_multislice.f90 \
	m_lens.f90 \
	m_tilt.f90 \
	quadpack.f90 \
	m_electron.f90 \
	m_absorption.f90 \
	m_potential.f90 \
	mod_Hn0.f90 \
	MS_utilities.f90

OTHER_SRC := \
	s_qep_tem.f90 \
	s_qep_stem.f90 \
	s_absorptive_tem.f90 \
	s_absorptive_stem.f90 \
	muSTEM.f90

ifeq ($(MODE),GPU)
  MOD_SRC := $(GPU_MOD_SRC)
else
  MOD_SRC := $(CPU_MOD_SRC)
endif

SRC := $(MOD_SRC) $(OTHER_SRC)
OBJ := $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))

# === Default Target ===
all: $(OUTDIR)/$(OUTNAME)

# === Compile Rule ===
$(OBJDIR)/%.o: %.f90
	@mkdir -p $(dir $@)
	@echo "Compiling $< → $@"
ifeq ($(MODE),GPU)
	$(FC) $(FCFLAGS) $(INCLUDES) -module $(OBJDIR) -c $< -o $@
else
	$(FC) $(FCFLAGS) $(INCLUDES) -module $(OBJDIR) -c $< -o $@
endif

# === Linking ===
$(OUTDIR)/$(OUTNAME): $(OBJ)
	@mkdir -p $(OUTDIR)
	@echo "Linking $@"
ifeq ($(MODE),GPU)
	$(FC) $(FCFLAGS) $(LIBDIRS) $(OBJ) $(LIBS) -o $@
else
	$(FC) $(FCFLAGS) $(LIBDIRS) $(OBJ) $(LIBS) -o $@
endif

# === Cleanup ===
clean:
	rm -rf $(OBJDIR) $(OUTDIR)/$(OUTNAME)

rebuild: clean all

.PHONY: all clean rebuild
