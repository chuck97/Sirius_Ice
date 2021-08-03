AFTHOME  := $(shell cd ../../;pwd)

all: help

lib:
	@if test ! -d freeglut-2.8.1 ; then \
	tar -xvzf freeglut-2.8.1.tar.gz; \
	cd freeglut-2.8.1; \
	./configure --prefix=$(AFTHOME) --disable-warnings --disable-dependency-tracking; \
	cd ..; \
	fi; \
	cd freeglut-2.8.1; \
	make; \
	make install;

clean:
	@if test -d freeglut-2.8.1 ; then \
	cd freeglut-2.8.1; \
	make clean; \
	fi;

purge: clean
	@if test -d freeglut-2.8.1 ; then \
	cd freeglut-2.8.1; \
	make uninstall; \
	cd ..; \
	rm -r freeglut-2.8.1; \
	fi;

help:
	@echo "make {lib|clean|purge|help}"
	@echo " "
	@echo "   lib   - installs the package freeglut"
	@echo "   clean - clean the package"
	@echo "   help  - print this message"
	@echo " "

