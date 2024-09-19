MWRAP_INSTALL = ${HOME}/mwrap/mwrap
FMM3DBIE_STATIC_INSTALL = ${HOME}/fmm3dbie/lib-static/libfmm3dbie_matlab.a
MAGNETOSTATICS = ${HOME}/fmm_taylor/virtual-casing-helmholtz/src/magneto-static-routs.o
MAGNETODYNAMICS = ${HOME}/fmm_taylor/virtual-casing-helmholtz/src/magneto-dynamic-routs.o
HELPER = ${HOME}/fmm_taylor/virtual-casing-helmholtz/src/fmm-helper-routs.o

# MWRAP_INSTALL = ${HOME}/Documents/fmm_taylor/mwrap/mwrap
# FMM3DBIE_STATIC_INSTALL = ${HOME}/Documents/fmm_taylor/fmm3dbie/lib-static/libfmm3dbie_matlab.a
# MAGNETOSTATICS = ${HOME}/Documents/fmm_taylor/virtual-casing/src/magneto-static-routs.o
# MAGNETODYNAMICS = ${HOME}/Documents/fmm_taylor/virtual-casing/src/magneto-dynamic-routs.o
# HELPER = ${HOME}/Documents/fmm_taylor/virtual-casing-helmholtz/src/fmm-helper-routs.o

matlab:
	mkdir -p +taylor/+static/
	$(MWRAP_INSTALL) -c99complex -list -mex gradcurlS0 -mb gradcurlS0.mw
	$(MWRAP_INSTALL) -c99complex -mex gradcurlS0 -c gradcurlS0.c gradcurlS0.mw
	mex gradcurlS0.c $(FMM3DBIE_STATIC_INSTALL) $(MAGNETOSTATICS) -compatibleArrayDims -DMWF77_UNDERSCORE1 "CFLAGS=-std=gnu17 -Wno-implicit-function-declaration -fPIC" -output gradcurlS0 -lm -lstdc++ -ldl -lgfortran -lgomp -lmwblas -lmwlapack
	mkdir -p +taylor/+dynamic/
	$(MWRAP_INSTALL) -c99complex -list -mex gradcurlSk -mb gradcurlSk.mw
	$(MWRAP_INSTALL) -c99complex -mex gradcurlSk -c gradcurlSk.c gradcurlSk.mw
	mex gradcurlSk.c $(FMM3DBIE_STATIC_INSTALL) $(MAGNETODYNAMICS) -compatibleArrayDims -DMWF77_UNDERSCORE1 "CFLAGS=-std=gnu17 -Wno-implicit-function-declaration -fPIC" -output gradcurlSk -lm -lstdc++ -ldl -lgfortran -lgomp -lmwblas -lmwlapack
	mkdir -p +taylor/+helper/
	$(MWRAP_INSTALL) -c99complex -list -mex helper -mb helper.mw
	$(MWRAP_INSTALL) -c99complex -mex helper -c helper.c helper.mw
	mex -v helper.c $(FMM3DBIE_STATIC_INSTALL) $(HELPER) $(MAGNETODYNAMICS) -compatibleArrayDims -DMWF77_UNDERSCORE1 "CFLAGS=-std=gnu17 -Wno-implicit-function-declaration -fPIC" -output helper -lm -lstdc++ -ldl -lgfortran -lgomp -lmwblas -lmwlapack

