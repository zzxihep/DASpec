detected_OS := $(shell uname)

ifeq ($(detected_OS), Linux)
	CFLAGS = -fPIC -O2 -Wall
	CC = g++
else
	CFLAGS = -fPIC -O2 -Wall -undefined dynamic_lookup
	CC = clang++
	# CC = g++-13
endif


suffix := $(shell python3-config --extension-suffix)
PY_CFLAGS := $(shell python3-config --includes)
# PYBIND11 := $(shell python3 -m pybind11 --includes)
GSL := $(shell pkg-config --cflags --libs gsl)
SHARE = -shared
MPFIT = -I ../mpfit/cmpfit-1.3a -L ../mpfit/cmpfit-1.3a -lmpfit
SWIG = swig -c++ -python


objects = component.o function.o compcontainer.o curvefit.o swigDASpec_wrap.o

default : _swigDASpec.so _carray.so

_swigDASpec.so : $(objects) common.h prototype.h DASpec.h
	$(CC) -o _swigDASpec.so $(objects) $(MPFIT) $(GSL) $(SHARE) $(PY_CFLAGS) $(CFLAGS)

component.o : component.cpp component.h prototype.h common.h function.h
	$(CC) -c component.cpp $(CFLAGS) $(GSL)

function.o : function.cpp function.h common.h
	$(CC) -c function.cpp $(CFLAGS) $(GSL)

compcontainer.o : compcontainer.cpp compcontainer.h common.h
	$(CC) -c compcontainer.cpp $(CFLAGS) $(GSL)

curvefit.o : curvefit.cpp curvefit.h common.h
	$(CC) -c curvefit.cpp $(CFLAGS) $(MPFIT) $(GSL)

swigDASpec_wrap.o : swigDASpec_wrap.cxx
	$(CC) -c swigDASpec_wrap.cxx $(CFLAGS) $(PY_CFLAGS) $(MPFIT)

swigDASpec_wrap.cxx : swigDASpec.i
	$(SWIG) swigDASpec.i


_carray.so : carray.o carray_wrap.o
	$(CC) -o _carray.so carray.o carray_wrap.o $(SHARE) $(PY_CFLAGS) $(CFLAGS)

carray.o : carray.cpp
	$(CC) -c carray.cpp $(CFLAGS)

carray_wrap.o : carray_wrap.cxx
	$(CC) -c carray_wrap.cxx $(CFLAGS) $(PY_CFLAGS)

carray_wrap.cxx : carray.i
	$(SWIG) carray.i


clean :
	rm _swigDASpec.so $(objects) swigDASpec_wrap.cxx carray_wrap.cxx carray.o carray_wrap.o _carray.so swigDASpec.py carray.py