############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make

LIBS    = -L$(ANILIB) \
	  $(LIBC2F) $(LIBMBA) $(LIBLMR) \
	  $(LIBLAPACK) $(LIBBLAS) -lm $(LIBSYS)


############################################################
EXEANALYTIC = $(ANIBIN)/mbaAnalytic.exe
EXENODAL    = $(ANIBIN)/mbaNodal.exe

OBJANALYTIC = main_analytic.o
OBJNODAL    = main_nodal.o


############################################################
all: 	    help

exe: 	    $(EXEANALYTIC) $(EXENODAL) info

run-ana:
	    cd $(ANIBIN) && ./mbaAnalytic.exe

run-nod:
	    cd $(ANIBIN) && ./mbaNodal.exe 

clean: 
	    @/bin/rm -f *.o  *~ $(ANIBIN)/*.gmv $(ANIBIN)/*.ani\
            $(ANIBIN)/*.exe $(ANIBIN)/gmon.out 

help:
	    @echo "make {exe|run-ana|run-nod|clean|help}"
	    @echo " "
	    @echo "   exe - create two executables"
	    @echo "   run-ana - run code with an analytic metric"	
	    @echo "   run-nod - run code with a nodal metric"
	    @echo "   clean   - clean this directory"
	    @echo " "

info:
	    @echo " "
	    @echo "Executable files are located in ../../../bin/"
	    @echo "Libraries are located in ../../../lib/"
	    @echo "GMV figures are located in ../../../bin/"
	    @echo " "

gmv:
	    @cd $(ANIBIN); $(VIEWER) -i save.gmv


############################################################
$(EXEANALYTIC): $(OBJANALYTIC) 
	@echo ''
	@echo 'Creating executable bin/mbaAnalytic.exe' 
	@$(CLINKER) $(LDFLAGS) -o $(EXEANALYTIC) $(OBJANALYTIC) $(LIBMBA) $(LIBVIEW) $(LIBS)
	@echo ''

$(EXENODAL): $(OBJNODAL) 
	@echo ''
	@echo 'Creating executable bin/mbaNodal.exe' 
	@$(CLINKER) $(LDFLAGS) -o $(EXENODAL) $(OBJNODAL) $(LIBMBA) $(LIBVIEW) $(LIBS)
	@echo ''


############################################################
main_analytic.o   : $(ANIINC)/ani3D.h
main_nodal.o : $(ANIINC)/ani3D.h


############################################################
.f.o:
	@echo $(F77)' -c ' $*.f 
	@$(F77) $(FFLAGS) -c $*.f -o $*.o -I$(ANIINC)
.c.o:
	@echo $(CC)'  -c ' $*.c 
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o -I$(ANIINC)

