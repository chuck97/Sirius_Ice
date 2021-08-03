############################################################
# Use first  command if you have the LAPACK library
# Use second command to compily necessary LAPACK routines 
############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make


############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
LIBS	    = -L$(ANILIB) $(LIBFEM) $(LIBMBA) $(LIBLMR) $(LIBILU) \
	                  $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm 

EXE         = $(ANIBIN)/Elasticity.exe

OBJ         = main.o forlibfem.o 


############################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) &&  ./Elasticity.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.gmv $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run|gmv|clean|help}"
	    @echo " "
	    @echo "   exe   - compile the code and link with libraries "
	    @echo "           libmba3D, libfem3D, liblmr3D, libview3D, libilu, liblapack and libblas"
	    @echo "   run   - run adaptive computation"
	    @echo "   gmv   - show GMV mesh"
	    @echo "   clean - clean the package"
	    @echo "   help  - print this message"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in   bin"
	    @echo "Library files are located in      lib"
	    @echo "Postscript figures are located in bin"
	    @echo " "

gmv:
	    @echo " "
	    @echo "GMV figures: mesh.gmv"
	    @echo " "
	    cd $(ANIBIN); $(VIEWER) -i mesh.gmv



############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/Elasticity.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS)
	@echo ''


############################################################
main.o      : main.f  
forlibfem.o : forlibfem.f


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

