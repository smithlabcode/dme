# Copyright (C) 2025 Andrew D Smith
# Author: Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
#

PROGS = dme2
CXX = g++ -std=c++17
CXXFLAGS = -Wall
DEBUGFLAGS = -g

INCLUDEDIRS = smithlab_cpp
INCLUDEARGS = $(addprefix -I ,$(INCLUDEDIRS))

ifdef DEBUG
CXXFLAGS += -g
else
CXXFLAGS += -O3 -DNDEBUG
endif

all: $(PROGS)

$(PROGS): $(addprefix smithlab_cpp/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< $(INCLUDEARGS)

dme2: dme2.cpp dme_tcm_workspace.o dme_zoops_workspace.o CTSet.o \
	ScoringMatrix.o Matrix.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(INCLUDEARGS)

clean:
	@-rm -f $(PROGS) *.o
.PHONY: clean
