############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make


############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
LIBS = -L$(ANILIB) $(LIBFEM) $(LIBVIEW) $(LIBMBA) $(LIBLU) $(LIBINB) $(LIBLAPACK) $(LIBBLAS) -lm

EXE  = $(ANIBIN)/StokesNavier.exe

OBJ  =  main.o forlibfem.o forlibinb.o


############################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) &&  ./StokesNavier.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps $(ANIBIN)/*.exe $(ANIBIN)/*.gmv

help:
	    @echo "make {exe|run|gmv|clean|help}"
	    @echo " "
	    @echo "   exe   - compile the code and link with libraries"
	    @echo "   run   - run adaptive computation"
	    @echo "   gmv   - draw the mesh and solution"
	    @echo "   clean - clean the package"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in   bin"
	    @echo " "

gmv:
	    @echo " "
	    @echo "GMV Files:"
	    @echo " "
 
	    @cd $(ANIBIN); $(VIEWER) -i velocity_x.gmv; $(VIEWER) -i pressure.gmv


############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/StokesNavier.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS) 
	@echo ''


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o
