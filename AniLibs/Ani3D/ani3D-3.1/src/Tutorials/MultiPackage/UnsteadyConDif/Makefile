############################################################
ANIHOME = ../../../..
include $(ANIHOME)/src/Rules.make

LIBS	= -L$(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBMBA) $(LIBVIEW)  \
	              $(LIBILU) $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
EXE         = $(ANIBIN)/UnsConDif.exe

OBJ         = main.o forlibfem.o 


############################################################
all: 	    help

exe: 	    $(EXE) info

run:
	    cd $(ANIBIN) &&  ./UnsConDif.exe

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.ps $(ANIBIN)/*.exe 

help:
	    @echo "make {exe|run|clean|help|gmv}"
	    @echo " "
	    @echo "   exe - compile the code and link with libraries "
	    @echo " "
	    @echo "   run - run adaptive computation"
	    @echo " "
	    @echo "   gmv    - show  mesh and solution initial and final"
	    @echo "   clean - clean the package"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in   bin"
	    @echo "GMV figures are located in bin"
	    @echo " "

gmv:
	    @echo " "
	    @echo " "
	    @echo "GMV figures:"
	    @echo "   1. solution at T=0"
	    @echo "   2. solution at T=0.5"
	    @echo " "
	    cd $(ANIBIN); $(VIEWER) -i solution_initial.gmv; $(VIEWER) -i solution_final.gmv



############################################################
$(EXE): $(OBJ) 
	@echo ''
	@echo 'Creating executable bin/UnsConDif.exe' 
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

