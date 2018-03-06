################################################################################
# Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
#
# (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on : Jun 07, 2012 by andrade
# Last update: Jun 07, 2012 by andrade
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
# User lib includes
###############################################################################

DIRS = EDGEGEN FMATCH KDTREE LINKERN UTIL

###############################################################################
# Lib and include definitions
###############################################################################

include ../Makefile.inc

###############################################################################
# Build Rules
###############################################################################

all: $(DIRS)

.PHONY: all strip clean doc docclean depclean $(DIRS)
.SUFFIXES: .cpp .o

$(DIRS):
	make -C $@

clean:
	@echo "--> Cleaning compiled..."

	for dir in $(DIRS); do \
		make -C $$dir clean; \
	done
