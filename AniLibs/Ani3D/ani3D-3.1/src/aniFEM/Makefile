ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

########################################################################
LIBOBJS = fem3Dtet.o fem3Dface.o fem3Dtri.o fem3Derr.o fem3Dsub.o \
	  opIDEN.o opGRAD.o opDIV.o opCURL.o opDUDX.o \
	  template.o bc.o algebra.o assemble.o \
	  forlibfem.o isys.o dof.o error.o


########################################################################
all:     help

lib:	 $(LIBFEM) info

clean:
	 @/bin/rm -f  *.o  *~

help:
	 @echo "make {lib|clean|help}"
	 @echo " "
	 @echo "   lib   - build library libfem3D-$(version).a"
	 @echo "   clean - clean this directory"
	 @echo "   help  - print help message"
	 @echo " "

info:
	 @echo " "
	 @echo "Libriry version: " $(version)
	 @echo "Libraries are located in ../../lib/"
	 @echo " "


########################################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC) -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


$(LIBFEM): $(LIBOBJS)
	@rm -rf $(LIBFEM)
	@ar -r $(LIBFEM) $(LIBOBJS)


########################################################################
fem3Dtet.o : fem3Dtet.fd
fem3Derr.o : fem3Dtet.fd
fem3Dtri.o : fem3Dtet.fd fem3Dtri.fd assemble.fd
fem3Dsub.o : fem3Dtet.fd

opIDEN.o   : fem3Dtet.fd
opGRAD.o   : fem3Dtet.fd
opDIV.o    : fem3Dtet.fd
opCURL.o   : fem3Dtet.fd
opDUDX.o   : fem3Dtet.fd

assemble.o : fem3Dtet.fd assemble.fd
template.o : fem3Dtet.fd assemble.fd 
forlibfem.o: fem3Dtet.fd assemble.fd

user_lib.o : fem3Dtet.fd assemble.fd

dof.o      : fem3Dtet.fd
isys.o     : assemble.fd

