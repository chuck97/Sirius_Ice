############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make
THISDIR := $(shell /bin/pwd)


############################################################
LIBS     = -L$(ANILIB) $(LIBMBA) $(LIBVIEW) \
	               $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
#  DO NOT CHANGE THE COMMANDS BELOW THIS LINE
############################################################
EXE	= $(ANIBIN)/mainView.exe


############################################################
OBJ	= main.o 

############################################################
all:	    help

exe:        $(EXE) info

run:        $(EXE)
	    @cd $(ANIBIN) && ./mainView.exe

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.gmv 

help:		
	    @echo "make {exe|run|gmv|clean|help}"
	    @echo " "
	    @echo "   exe - create executable"
	    @echo "   run - load AFT mesh and save it in GMV format"
	    @echo "   gmv - draw the mesh and solution"
	    @echo "   clean   - clean the package"
	    @echo " "

info:	
	    @echo " " 
	    @echo "Summary:"
	    @echo "Executables are located in ../../../bin/"
	    @echo "Libraries  are  located in ../../../lib/"
	    @echo " "

gmv:	    
	    @cd $(ANIBIN); $(VIEWER) -a $(THISDIR)/gmv.attr -i save.gmv


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


$(EXE): $(OBJ)
	@echo ''
	@echo 'Creating executable $(EXE)'
	@$(FLINKER) $(LDFLAGS) -o $(EXE) $(OBJ) $(LIBS)
	@echo ''

