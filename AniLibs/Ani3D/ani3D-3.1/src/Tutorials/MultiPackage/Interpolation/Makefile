############################################################
# Use first  command if you have the LAPACK library
# Use second command to compily necessary LAPACK routines 
############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make


############################################################
# YOU DO NOT NEED TO CHANGE THE LINES BELOW
############################################################
LIBS	    = -L$(ANILIB) $(LIBMBA) $(LIBLMR) $(LIBVIEW) \
	                  $(LIBLAPACK) $(LIBBLAS) -lm

EXE         = $(ANIBIN)/Interpolation.exe

OBJ         = main.o esterr.o


############################################################
all: 	    help

exe: 	    $(EXE)

run:
	    cd $(ANIBIN) &&  ./Interpolation.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.gmv $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run|gmv|clean|help}"
	    @echo " "
	    @echo "     exe   - compile the code and link it with libraries"
	    @echo "             libmba3D, liblmr3D, libview3D, liblapack and libblas"
	    @echo "     run   - run adaptive computation"
	    @echo "     gmv   - show adaptive mesh save.gmv"
	    @echo "     clean - clean the package"
	    @echo " "

gmv:
	    @echo " "
	    @echo "GMV Viewer:"
	    @echo " "
	    cd $(ANIBIN); $(VIEWER) -i save.gmv


############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/Interpolation.exe' 
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS)
	@echo ''


############################################################
main.o      : main.f  
esterr.o    : esterr.f  


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o 
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o

