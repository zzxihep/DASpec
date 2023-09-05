detected_OS := $(shell uname)

ifeq ($(detected_OS), Linux)
	CFLAGS = -fPIC -O2 -Wall
	CC = g++
else
	CFLAGS = -fPIC -O2 -Wall -undefined dynamic_lookup
	# CC = clang++
	CC = g++-13
endif


suffix := $(shell python3-config --extension-suffix)
PY_CFLAGS := $(shell python3-config --includes)
# PYBIND11 := $(shell python3 -m pybind11 --includes)
GSL := $(shell pkg-config --cflags --libs gsl)
SHARE = -shared
MPFIT = -I ../mpfit/cmpfit-1.3a -L ../mpfit/cmpfit-1.3a -lmpfit


objects = component.o function.o compcontainer.o curvefit.o swigDASpec_wrap.o

default : _swigDASpec.so

_swigDASpec.so : $(objects) common.h prototype.h DASpec.h
	$(CC) -o _swigDASpec.so $(objects) $(MPFIT) $(GSL) $(SHARE) $(PY_CFLAGS) $(CFLAGS)

component.o : component.cpp component.h prototype.h common.h function.h
	$(CC) -c component.cpp $(CFLAGS)

function.o : function.cpp function.h common.h
	$(CC) -c function.cpp $(CFLAGS)

compcontainer.o : compcontainer.cpp compcontainer.h common.h
	$(CC) -c compcontainer.cpp $(CFLAGS)

curvefit.o : curvefit.cpp curvefit.h common.h
	$(CC) -c curvefit.cpp $(CFLAGS) $(MPFIT)

swigDASpec_wrap.o : swigDASpec_wrap.cxx
	$(CC) -c swigDASpec_wrap.cxx $(CFLAGS) $(PY_CFLAGS) $(MPFIT)

