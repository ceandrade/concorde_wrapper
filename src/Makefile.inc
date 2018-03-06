################################################################################
# Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
#
# (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jun 07, 2012 by andrade
# Last update: Jul 17, 2012 by andrade
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
# User Flags
###############################################################################

# Compiler flags for debugging
#USER_FLAGS = -g3 -Wall -Wextra -pedantic -fexceptions -fno-omit-frame-pointer \
#    -fno-optimize-sibling-calls -fno-inline

USER_FLAGS += -Wall -Wextra -pedantic -O3 -pthread -fopenmp \
    -funroll-loops -fpeel-loops -fprefetch-loop-arrays \
    -fomit-frame-pointer -ftracer

#LDFLAGS = -fdata-sections -ffunction-sections

###############################################################################
# Common libs
###############################################################################

LIBS = -lm -pthread -fopenmp

###############################################################################
# Compiler and linker defs
###############################################################################

# C Compiler command and flags
CC = gcc
CFLAGS = $(USER_FLAGS)

# C++ Compiler command and flags
CXX = g++
CXXFLAGS = $(USER_FLAGS)

# Lib maker commands
AR = ar
ARFLAGS	= rv
RANLIB = ranlib

# Other includes
RM = rm
SHELL = /bin/bash

