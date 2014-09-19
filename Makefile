#
# Makefile for non-Microsoft compilers
#

## Linux  (uncomment the 2 lines below for compilation on Linux)
CXXFLAGS += -Wall -std=c++11
LDFLAGS += -lrt -pthread -std=c++11

SNAPDIR = /home/jzzhao/git_project/netsnap
LIB = -I $(SNAPDIR)/glib -I $(SNAPDIR)/snap
OBJ = Snap.o TCEM.o

## Main application file
MAIN = main

all: $(MAIN)

opt: CXXFLAGS += -O4
opt: LDFLAGS += -O4
opt: $(MAIN)

# COMPILE
$(MAIN): $(MAIN).cpp stdafx.h $(OBJ) 
	g++ $(DEBUG) $(OBJ) $(LDFLAGS) $(LIB) -o $(MAIN) $<

TCEM.o: TCEM.cpp TCEM.h stdafx.h Snap.o
	g++ -c $(DEBUG) $(CXXFLAGS) $(LIB) $<

Snap.o: 
	g++ -c $(CXXFLAGS) $(LIB) $(SNAPDIR)/snap/Snap.cpp

clean:
	rm -f *.o $(MAIN) $(MAIN).exe $(MAIN).Err
	rm -rf Debug Release

