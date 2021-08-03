ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
OBJCORE     = ani3D.o


############################################################
all: 	    help

lib:	    $(LIBC2F) info

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.gmv $(ANIBIN)/*.exe  

help:
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "  lib   - build library libc2f3D-$(version).a"
	    @echo "  clean - clean this directory"
	    @echo "  help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Libriry version: " $(version)
	    @echo "Libraries are located in ../../lib"
	    @echo " "


############################################################
$(LIBC2F): $(OBJCORE) 
	@echo ''
	@rm -f $(LIBC2F)
	@ar -r $(LIBC2F) $(OBJCORE)


############################################################
ani3D.o         : ani3D.h


############################################################
.f.o:
	@echo $(F77)' -c ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo $(CC)'  -c ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

