############################################################
# Use first  command if you have the LAPACK library
# Use second command to compily necessary LAPACK routines 
############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBFEM) $(LIBMBA) $(LIBLMR) $(LIBILU) \
	              $(LIBLAPACK) $(LIBBLAS) $(LIBVIEW) -lm 


############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
EXE         = $(ANIBIN)/AnisotropicAdap.exe

OBJ         = main.o forlibfem.o 


############################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) &&  ./AnisotropicAdap.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.gmv $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run|gmv|clean|help}"
	    @echo " "
	    @echo "   exe   - compile the code and link with libraries "
	    @echo "           libmba3D, libfem3Di, liblmr3D, libview3D, libilu, liblapack and libblas"
	    @echo "   run   - run adaptive computation"
	    @echo "   gmv   - view GMV mesh"
	    @echo "   clean - clean the package"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in bin"
	    @echo "Library files are located in lib"
	    @echo "GMV files are located in bin"
	    @echo " "

gmv:
	    @echo " "
	    @echo "GMV figures:"
	    @echo " "
	    cd $(ANIBIN); $(VIEWER) -i mesh.gmv



############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/AnisotropicAdap.exe' 
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

