matlab:
	${HOME}/mwrap/mwrap -c99complex -list -mex gradcurlS0 -mb gradcurlS0.mw
	${HOME}/mwrap/mwrap -c99complex -mex gradcurlS0 -c gradcurlS0.c gradcurlS0.mw
	mex gradcurlS0.c ${HOME}/fmm3dbie/lib-static/libfmm3dbie.a ${HOME}/fmm_taylor/virtual-casing-helmholtz/src/magneto-static-routs.o -compatibleArrayDims -DMWF77_UNDERSCORE1 "CFLAGS=-std=gnu17 -Wno-implicit-function-declaration -fPIC" -output gradcurlS0 -lm -lstdc++ -ldl -lgfortran -lgomp -lmwblas -lmwlapack
