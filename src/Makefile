########################################################################################
#
#  KPP - The Kinetic PreProcessor
#        Builds simulation code for chemical kinetic systems
#
#  Copyright (C) 1995-1996 Valeriu Damian and Adrian Sandu
#  Copyright (C) 1997-2005 Adrian Sandu
#      with contributions from: Mirela Damian, Rolf Sander 
#
#  KPP is free software; you can redistribute it and/or modify it under the
#  terms of the GNU General Public License as published by the Free Software
#  Foundation (http://www.gnu.org/copyleft/gpl.html); either version 2 of the
#  License, or (at your option) any later version.
#
#  KPP is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#  details.
#
#  You should have received a copy of the GNU General Public License along
##  with this program; if not, consult http://www.gnu.org/copyleft/gpl.html or
#  write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
#  Boston, MA  02111-1307,  USA.
#
#  Adrian Sandu
#  Computer Science Department
#  Virginia Polytechnic Institute and State University
#  Blacksburg, VA 24060
#  E-mail: sandu@cs.vt.edu
#
#######################################################################################

# Load machine-specific settings for CC, CC_FLAGS etc.:
include Makefile.defs

PROG = ../bin/kpp

OBJS  = \
	code.o \
	code_c.o \
	code_f77.o \
	code_f90.o \
	code_matlab.o \
	debug.o \
	gen.o \
	kpp.o \
	lex.yy.o \
	scanner.o \
	scanutil.o \
	y.tab.o

.PHONY: all
all: $(OBJS)
	$(CC) $(CC_FLAGS) $(OBJS) -L$(FLEX_LIB_DIR) -lfl -o kpp
	@mkdir -p ../bin
	@mv kpp $(PROG)

.PHONY: clean
clean:  
	@rm -f *~ *.o 

.PHONY: distclean
distclean: clean
	@rm -f $(PROG)

.PHONY: maintainer-clean
maintainer-clean: distclean
	@rm -f lex.yy.c y.tab.c y.tab.h

# list the configuration:
.PHONY: list
list:
	@echo "SYSTEM       = $(SYSTEM)"
	@echo "HOST         = $(HOST)"
	@echo "CC           = $(CC)"
	@echo "CC_FLAGS     = $(CC_FLAGS)"
	@echo "FLEX         = $(FLEX)"
	@echo "FLEX_LIB_DIR = $(FLEX_LIB_DIR)"
	@echo "BISON        = $(BISON)"
	@echo "INCLUDE_DIR  = $(INCLUDE_DIR)"

.c.o:
	$(CC) $(CC_FLAGS) -I$(INCLUDE_DIR) -c $*.c

lex.yy.c: scan.l scan.h gdata.h
	$(FLEX) -olex.yy.c scan.l 

y.tab.c: scan.y scan.h
	$(BISON) -d -o y.tab.c scan.y

y.tab.h: scan.y scan.h
	$(BISON) -d -o y.tab.c scan.y

# dependencies via header files:
code.o:        code.h gdata.h
code_c.o:      code.h gdata.h
code_f77.o:    code.h gdata.h
code_f90.o:    code.h gdata.h
code_matlab.o: code.h gdata.h
debug.o:              gdata.h scan.h
gen.o:         code.h gdata.h scan.h
kpp.o:                gdata.h scan.h
lex.yy.o:             gdata.h scan.h y.tab.h
scanner.o:            gdata.h scan.h 
scanutil.o:           gdata.h scan.h
y.tab.o:                      scan.h y.tab.h
