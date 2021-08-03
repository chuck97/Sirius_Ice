############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
ANIHOME = ../..
include $(ANIHOME)/src/Rules.make


############################################################
OBJVIEW =  drawGMV.o


############################################################
all: 	    help

lib:	    $(LIBVIEW) info

clean: 
	    @/bin/rm -f *.o  

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib   - create library libview3D-$(version).a"
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Libriry version: " $(version)
	    @echo "Libraries are located in ../../lib"
	    @echo " "


############################################################
$(LIBVIEW): $(OBJVIEW)
	@echo ''
	@rm -f $(LIBVIEW)
	@ar -r $(LIBVIEW) $(OBJVIEW)


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

