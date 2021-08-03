ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
LIBOBJS  = comp_geom.o  metamesh_fast.o  metamesh_slow.o rhs.o  tet2tet.o

############################################################
all:	    help

lib:	    $(LIBPRJ) info

clean:
	    @/bin/rm -f *.o *~ 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - build library libprj3D-$(version).a"
	    @echo " "
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo "Summary:"
	    @echo "Libraries   are located in lib/"
	    @echo " "

############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)


############################################################
$(LIBPRJ): $(LIBOBJS)
	   @echo ''
	   @rm -f $(LIBPRJ)
	   @ar -r $(LIBPRJ)  $(LIBOBJS)

############################################################

# rhs.o : fem3Dtet.fd assemble.fd

