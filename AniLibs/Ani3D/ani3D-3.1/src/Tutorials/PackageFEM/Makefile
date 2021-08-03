############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS    = -L$(ANILIB) $(LIBFEM) $(LIBMBA) $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
EXEFEMTET = $(ANIBIN)/aniFEM_tetrahedron.exe 
EXEFEMBC  = $(ANIBIN)/aniFEM_bc.exe
EXETEMFEM = $(ANIBIN)/aniFEM_template.exe

FEMTETOBJS = main_tetrahedron.o
FEMOBJS    = main_simple.o  
FEMBCOBJS  = main_bc.o  
FEMTEMOBJS = main_template.o


############################################################
all:     help

exe:     $(EXEFEMTET) $(EXEFEM) $(EXETEMFEM) $(EXEFEMBC) info

run-tet:
	 @cd $(ANIBIN) && ./aniFEM_tetrahedron.exe

run-bc:  $(EXEFEMBC)
	 @cd $(ANIBIN) && ./aniFEM_bc.exe

run-tem: $(EXETEMFEM)
	 @cd $(ANIBIN) && ./aniFEM_template.exe

clean:
	 @/bin/rm -f  *.o  *~  $(ANIBIN)/*.exe $(ANIBIN)/CSRsystem $(ANIBIN)/*.gmv

help:
	 @echo "make {exe|run-bc|run-tem|clean|help}"
	 @echo " "
	 @echo "   exe     - compile three executabels"
	 @echo "   run-tet - generation of elemental matrices" 
	 @echo "   run-tem - run code with templated local matrix"
	 @echo "   run-bc  - generate matrix for diffusion problem with various boundary conditions"
	 @echo "   clean   - clean this directory"
	 @echo " "

info:
	 @echo " "
	 @echo "Executable files are located in ../../../bin/"
	 @echo " "


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC) -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


############################################################
$(EXEFEMTET): $(FEMTETOBJS)
	@echo 'Creating executable $(EXEFEMTET)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEFEMTET) $(FEMTETOBJS) $(LIBS)

$(EXEFEMBC): $(FEMBCOBJS)
	@echo 'Creating executable $(EXEFEMBC)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEFEMBC) $(FEMBCOBJS) $(LIBS)

$(EXETEMFEM): $(FEMTEMOBJS)
	@echo 'Creating executable $(EXETEMFEM)'
	@$(FLINKER) $(LDFLAGS) -o $(EXETEMFEM) $(FEMTEMOBJS) $(LIBS)


############################################################



