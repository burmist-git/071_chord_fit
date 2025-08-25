ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX  = g++
CXX += -I./

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated
CXXFLAGS += $(ROOTCFLAGS)
CXXFLAGS += $(ROOTLIBS)
CXXFLAGS += $(ROOTGLIBS)
CXXFLAGS += -fconcepts
CXXFLAGS += -lMinuit

#----------------------------------------------------#

all: D_cord

.PHONY: printmakehelp_and_reminder
printmakehelp_and_reminder: D_cord.cpp Makefile
	$(info  /***********************************************************/)
	$(info  * task --> printmakehelp_and_reminder: D_cord.cpp Makefile *)
	$(info  * $$@ ----> $@                      *)
	$(info  * $$< --------------------------------> $<          *)
	$(info  * $$^ --------------------------------> $^ *)
	$(info  /***********************************************************/)

D_cord: D_cord.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm -rf D_cord
	rm -rf *~
