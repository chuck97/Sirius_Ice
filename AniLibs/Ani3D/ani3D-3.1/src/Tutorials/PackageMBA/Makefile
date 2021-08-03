############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS     = -L$(ANILIB) $(LIBMBA) $(LIBLMR) \
	               $(LIBLAPACK) $(LIBBLAS) -lm


############################################################
#  DO NOT CHANGE THE COMMANDS BELOW THIS LINE
############################################################
EXEANA	= $(ANIBIN)/mbaAnalytic.exe
EXENOD	= $(ANIBIN)/mbaNodal.exe

EXEAFT  = $(ANIBIN)/mbaReadAft.exe
EXEFIX  = $(ANIBIN)/mbaFixAft.exe
EXEREF  = $(ANIBIN)/mbaRefine.exe
EXEXTR  = $(ANIBIN)/mbaExtractFixUnite.exe


############################################################
OBJANA	= main_analytic.o 
OBJNOD	= main_nodal.o 

OBJAFT  = main_read_aft_mesh.o
OBJFIX  = main_fixshape.o
OBJREF  = main_refine.o
OBJXTR  = main_split_unite.o

############################################################
all:	    help

exe:        $(EXEANA) $(EXENOD) $(EXEAFT) $(EXEFIX) $(EXEREF) $(EXEXTR) info

run-ana:    $(EXEANA)
	    @cd $(ANIBIN) && ./mbaAnalytic.exe

run-nod:    $(EXENOD)
	    @cd $(ANIBIN) && ./mbaNodal.exe

run-aft:    $(EXEAFT)
	     @cd $(ANIBIN) && ./mbaReadAft.exe

run-fix:    $(EXEAFT)
	     @cd $(ANIBIN) && ./mbaFixAft.exe

run-ref:    $(EXEREF)
	     @cd $(ANIBIN) && ./mbaRefine.exe

run-xtr:    $(EXEXTR)
	     @cd $(ANIBIN) && ./mbaExtractFixUnite.exe

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.gmv $(ANIBIN)/*.nodal $(ANIBIN)/*.ani $(ANIBIN)/save.*

help:		
	    @echo "make {exe|run-ana|run-nod|run-aft|run-fix|run-ref|gmv|clean|help}"
	    @echo " "
	    @echo "   exe     - create five executables"
	    @echo " "
	    @echo "   run-ana - run code with analytic user-defined metric"
	    @echo "   run-nod - run code with nodal user-defined function"
	    @echo "   run-aft - load AFT mesh (.out) and save it as MBA (.ani) mesh"
	    @echo "   run-fix - load AFT mesh, fix its shape and save it as AFT mesh"
	    @echo "   run-ref - load and refine the mesh"
	    @echo "   run-xtr - load mesh, extract submeshes, fix their shape, unite submeshes, save  as AFT mesh"
	    @echo " "
	    @echo "   gmv     - visualize the final mesh"
	    @echo "   clean   - clean the package"
	    @echo " "

info:	
	    @echo " " 
	    @echo "Summary:"
	    @echo "Executables are located in ../../../bin/"
	    @echo "Libraries  are  located in ../../../lib/"
	    @echo " "

gmv:	    
	    @cd $(ANIBIN); $(VIEWER) -i save.gmv


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


$(EXEANA): $(OBJANA)
	@echo ''
	@echo 'Creating executable $(EXEANA)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEANA) $(OBJANA) $(LIBS)
	@echo ''

$(EXENOD): $(OBJNOD)
	@echo ''
	@echo 'Creating executable $(EXENOD)'
	@$(FLINKER) $(LDFLAGS) -o $(EXENOD) $(OBJNOD) $(LIBS)
	@echo ''

$(EXEAFT): $(OBJAFT)
	@echo ''
	@echo 'Creating executable $(EXEAFT)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEAFT) $(OBJAFT) $(LIBS)
	@echo ''

$(EXEFIX): $(OBJFIX)
	@echo ''
	@echo 'Creating executable $(EXEFIX)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEFIX) $(OBJFIX) $(LIBS)
	@echo ''

$(EXEREF): $(OBJREF)
	@echo ''
	@echo 'Creating executable $(EXEREF)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEREF) $(OBJREF) $(LIBS)
	@echo ''

$(EXEXTR): $(OBJXTR)
	@echo ''
	@echo 'Creating executable $(EXEXTR)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEXTR) $(OBJXTR) $(LIBS)
	@echo ''


