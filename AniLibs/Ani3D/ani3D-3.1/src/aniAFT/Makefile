ANIHOME = ../..
#SYSHOME := $(shell cd ../..; /bin/pwd)
include $(ANIHOME)/src/Rules.make


############################################################
AFTSRC  = src


############################################################
all: 	
	@echo "make {lib|gui|clean|purge}"
	@echo " "
	@echo "   lib   - build libraries in ../../lib"
	@echo "   gui   - create glut-based ../../bin/gui_glut.exe"
	@echo "   clean - clean"
	@echo "   purge - clean + remove binaries and libs"

gui: lib
	+cd src/aniFRT/GUI && make exe

lib:
	+cd src && make lib

clean: 
	+cd src && make clean

purge:
	+cd src && make purge
