#
# Makefile for non-Microsoft compilers
#

## Linux  (uncomment the 2 lines below for compilation on Linux)
CXXFLAGS += -Wall -std=c++11
LDFLAGS += -lrt -pthread -std=c++11

SNAPDIR = /home/jzzhao/git_project/netsnap
LIB = -I $(SNAPDIR)/glib -I $(SNAPDIR)/snap
OBJ = Snap.o TCEM.o ExamMgr.o TCEMGeneral.o

all: th

opt: CXXFLAGS += -O4
opt: LDFLAGS += -O4

# COMPILE
th: main_th.cpp $(OBJ) 
	g++ $(DEBUG) $(OBJ) $(LDFLAGS) $(LIB) -o $@ $<

nth: main_nth.cpp $(OBJ) 
	g++ $(DEBUG) $(OBJ) $(LDFLAGS) $(LIB) -o $@ $<
	
st: main_st.cpp $(OBJ) 
	g++ $(DEBUG) $(OBJ) $(LDFLAGS) $(LIB) -o $@ $<

test: main_test.cpp $(OBJ) 
	g++ $(DEBUG) $(OBJ) $(LDFLAGS) $(LIB) -o $@ $<
	
TCEMGeneral.o: TCEMGeneral.cpp TCEMGeneral.h stdafx.h Snap.o
	g++ -c $(DEBUG) $(CXXFLAGS) $(LIB) $<
	
TCEM.o: TCEM.cpp TCEM.h stdafx.h Snap.o
	g++ -c $(DEBUG) $(CXXFLAGS) $(LIB) $<

ExamMgr.o: ExamMgr.cpp ExamMgr.h stdafx.h Snap.o
	g++ -c $(DEBUG) $(CXXFLAGS) $(LIB) $<

Snap.o: $(SNAPDIR)/snap/Snap.cpp
	g++ -c $(CXXFLAGS) $(LIB) $<

clean:
	rm -f *.o th nth st test *.Err
	rm -rf Debug Release

