#
# Makefile for non-Microsoft compilers
#

## Linux  (uncomment the 2 lines below for compilation on Linux)
CXXFLAGS += -Wall -std=c++11
LDFLAGS += -lrt -pthread -std=c++11

SNAPDIR = /home/jzzhao/git_project/netsnap
LIB = -I $(SNAPDIR)/glib -I $(SNAPDIR)/snap
OBJ = Snap.o TCEM.o ExamMgr.o TCEMGeneral.o

all: main_th
th: main_th
nth: main_nth
st: main_st

opt: CXXFLAGS += -O4
opt: LDFLAGS += -O4

# COMPILE
main_th: main_th.cpp stdafx.h $(OBJ) 
	g++ $(DEBUG) $(OBJ) $(LDFLAGS) $(LIB) -o $@ $<

main_nth: main_nth.cpp stdafx.h $(OBJ) 
	g++ $(DEBUG) $(OBJ) $(LDFLAGS) $(LIB) -o $@ $<
	
main_st: main_st.cpp stdafx.h $(OBJ) 
	g++ $(DEBUG) $(OBJ) $(LDFLAGS) $(LIB) -o $@ $<

TCEMGeneral.o: TCEMGeneral.cpp TCEMGeneral.h stdafx.h Snap.o
	g++ -c $(DEBUG) $(CXXFLAGS) $(LIB) $<
	
TCEM.o: TCEM.cpp TCEM.h stdafx.h Snap.o
	g++ -c $(DEBUG) $(CXXFLAGS) $(LIB) $<

ExamMgr.o: ExamMgr.cpp ExamMgr.h stdafx.h Snap.o
	g++ -c $(DEBUG) $(CXXFLAGS) $(LIB) $<

Snap.o: 
	g++ -c $(CXXFLAGS) $(LIB) $(SNAPDIR)/snap/Snap.cpp

clean:
	rm -f *.o $(MAIN) $(MAIN).exe $(MAIN).Err
	rm -rf Debug Release

