################################################################################
# Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
#
# (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jul 09, 2011 by andrade
# Last update: Jul 31, 2012 by andrade
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################

###############################################################################
# User Defines
###############################################################################

# Set debug mode (you must set flags on Makefile.inc too).
#USER_DEFINES = -DDEBUG -DFULLDEBUG

include Makefile.inc

###############################################################################
# Build options
###############################################################################

EXE = test

###############################################################################
# The user source files
###############################################################################

# Include dir
USER_INCDIRS = . \
	./concorde/INCLUDE \
	./concorde_wrapper

# Source directories (to clean up object files)
SRC_DIRS = . \
	./concorde_wrapper

# Object files
OBJS = \
	./concorde_wrapper/concorde_wrapper.o

# Main function
OBJS += ./test_tsplib_loading.o

###############################################################################
# Lib and include definitions
###############################################################################

##############################
# Concorde libs
##############################

CONCORDE_PATH = ./concorde
CONCORDE_LIBS = \
	$(CONCORDE_PATH)/LINKERN/linkern.a \
	$(CONCORDE_PATH)/EDGEGEN/edgegen.a \
	$(CONCORDE_PATH)/KDTREE/kdtree.a \
	$(CONCORDE_PATH)/FMATCH/fmatch.a \
	$(CONCORDE_PATH)/UTIL/util.a

##############################
# Consolidate paths
##############################

# Consolidate include paths
USER_INCDIRS += /opt/local/include
	
# Consolidate include paths
USER_LIBDIRS += /opt/local/lib

# Libraries necessary to link.
LIBS += -lemon
	
###############################################################################
# Compiler defs
###############################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = $(USER_FLAGS)

# Necessary Include dirs
# Put -I in front of dirs
INCLUDES = `for i in $(USER_INCDIRS); do echo $$i | sed -e 's/^/-I/'; done`

# Necessary library dirs
# Put -L in front of dirs
LIBDIRS = `for i in $(USER_LIBDIRS); do echo $$i | sed -e 's/^/-L/'; done`

###############################################################################
# Build Rules
###############################################################################

all: $(EXE)

.PHONY: all concorde strip clean doc docclean depclean
.SUFFIXES: .cpp .o

$(EXE): concorde $(OBJS)
	@echo "--> Linking objects... "
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(CONCORDE_LIBS) $(MAURO_LIB) $(LIBDIRS) $(LIBS) -o $@
	@echo

.cpp.o:
	@echo "--> Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(USER_DEFINES) -c $< -o $@
	@echo
	
concorde:
	@echo "|---------------------[ Compiling Concorde ]---------------------|"
	make -C concorde
	@echo "|-----------------[ End of compiling Concorde ]------------------|"

strip: $(EXE)
	@echo "-->Stripping $<..."
	strip $<
	@echo

doc:
	doxygen Doxyfile

clean:
	@echo "--> Cleaning compiled..."
	rm -f $(EXE) $(OBJS)
	rm -rf Debug
	make -C concorde clean
	
depclean: clean docclean
	rm -rf `for i in $(SRC_DIRS); do echo $$i*~; done` 
		
docclean:
	@echo "--> Cleaning doc..."
	rm -rf doc

