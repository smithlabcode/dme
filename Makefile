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
LIBS = -lpopt
CXX = g++
CFLAGS = -Wall
DEBUGFLAGS = -g
INCLUDEDIR = cread
LIBDIR = cread

CREAD_SOURCES = $(shell ls cread/*.cpp)
CREAD_OBJS = $(subst .cpp,.o,$(CREAD_SOURCES))
OPT = 1

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CFLAGS += -O2
endif

all:	$(PROGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ -I$(INCLUDEDIR) $<

dme2:	dme2.o dme_tcm_workspace.o dme_zoops_workspace.o CTSet.o
	make -C cread	
	$(CXX) $(CFLAGS) -o $@ $^ $(CREAD_OBJS) $(LIBS) -I$(INCLUDEDIR) -L$(LIBDIR)

clean:
	make -C cread clean
	@-rm -f $(PROGS) *.o
.PHONY: clean
