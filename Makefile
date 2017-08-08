#
# Copyright (C) 2008 Cold Spring Harbor Laboratory and Andrew D Smith
# Author: Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# 

PROGS = dme2
CXX = g++ -std=c++11
CFLAGS = -Wall -O2
DEBUGFLAGS = -g

# check if a global copy of smithlab_cpp cannot be found and try to
# use a copy that is in the current directory
ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(abspath $(dir $(MAKEFILE_LIST)))/smithlab_cpp
ifeq ("$(wildcard $(SMITHLAB_CPP))","")
$(error SMITHLAB_CPP not set and smithlab_cpp not found)
endif
endif

INCLUDEDIRS = $(SMITHLAB_CPP)
INCLUDEARGS = $(addprefix -I ,$(INCLUDEDIRS))

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

all:	$(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

dme2:	dme2.cpp dme_tcm_workspace.o dme_zoops_workspace.o CTSet.o \
	Pattern.o Motif.o ScoringMatrix.o MotifSite.o Matrix.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS) $(INCLUDEARGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $(INCLUDEARGS) $<

clean:
	@-rm -f $(PROGS) *.o
.PHONY: clean
