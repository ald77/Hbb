CXX = $(shell root-config --cxx)
CXXFLAGS = -isystem$(shell root-config --incdir) -Wall -Wextra -pedantic -Wshadow $(shell root-config --cflags) -O2
LD = $(shell root-config --ld)
LDFLAGS = $(shell root-config --ldflags)
LDLIBS = $(shell root-config --libs) -lMinuit

default: n_minus_1.exe

n_minus_1.exe: n_minus_1.o cut.o
	$(LD) $(LDFLAGS) -o n_minus_1.exe n_minus_1.o cut.o $(LDLIBS)

n_minus_1.o: n_minus_1.cpp cut.hpp my_style.hpp
	$(CXX) $(CXXFLAGS) -c n_minus_1.cpp


.PHONY: clean

clean:
	-rm -rf *.exe *.o *~ *#
