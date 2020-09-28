# Makefile for the ROOT test progrrpc.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linuxegcs

CXX           =
ObjSuf        = o
SrcSuf        = cpp
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 

EVENTLIB      = $(EVENTSO)

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


ifeq ($(ARCH),linuxegcs)
CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC -g
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------
EXEFILO       = ProcessSignal.$(ObjSuf)
EXEFILES      = ProcessSignal.$(SrcSuf)
EXEFILE       = getParas$(ExeSuf)
GETFILTERTO   = GetFilter.$(ObjSuf)
GETFILTERTS   = GetFilter.$(SrcSub)

OBJS          = $(EXEFILO) $(GETFILTERTO)
PROGRRPC      = $(EXEFILE)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRRPC)

$(EXEFILE):     $(EXEFILO) $(GETFILTERTO)
		$(LD) $(LDFLAGS) $(EXEFILO) $(GETFILTERTO)\
		$(GLIBS) -lMinuit $(OutPutOpt)$(EXEFILE)
		@echo "$@ done"

clean:
		@rm -f $(OBJS) core $(EXEFILE) *~

distclean:      clean
		@rm -f $(PROGRRPC) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)

###
$(GETFILTERTO): GetFilter.h GetFilter.cc
$(EXEFILO): ProcessSignal.cc GetFilter.cc GetFilter.h
.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
