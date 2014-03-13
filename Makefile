# For debugging purposes, set -O0 and -g.
# For production use, set -O3.
CXXFLAGS = -O3 -g -Wall -pedantic -lOpenCL -fopenmp
LDLIBS = 
CXX = g++
COMPILE.cc = ${CXX} ${CXXFLAGS} ${CPPFLAGS} -c

SOURCES := $(wildcard [^_]*.cc)
OBJECTS := ${SOURCES:.cc=.o}
BINARY = raytrace

all: $(OBJECTS)
	$(CXX) -lOpenCL -fopenmp -msse -msse2 -msse3 -o $(BINARY) $(OBJECTS) $(LDLIBS)

%.o: %.cc
	${COMPILE.cc} -o $@ $<

clean:
	$(RM) $(BINARY) $(OBJECTS)
