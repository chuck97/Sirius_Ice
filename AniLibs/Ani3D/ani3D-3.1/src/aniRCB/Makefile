ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
LIBOBJS  = aux.o refine.o coarse.o 

############################################################
all:	    help

lib:	    $(LIBRCB) info

clean:
	    @/bin/rm -f *.o *~ 

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib - build library librcb3D-$(version).a"
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
	@$(F77) $(FFLAGS) -c $*.f -o $*.o 


############################################################
$(LIBRCB): $(LIBOBJS)
	   @echo ''
	   @rm -f $(LIBRCB)
	   @ar -r $(LIBRCB)  $(LIBOBJS)


############################################################
