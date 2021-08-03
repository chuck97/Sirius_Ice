############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make
THISDIR := $(shell /bin/pwd)

############################################################
LIBS     = -L$(ANILIB) $(LIBPRJ) $(LIBFEM) $(LIBMBA) $(LIBILU) \
           $(LIBLAPACK) $(LIBBLAS) $(LIBVIEW) -lm 

############################################################
EXEPRJ	= $(ANIBIN)/aniPRJ.exe
OBJPRJ	= main.o  

############################################################
all:	    help

exe:        $(EXEPRJ) info

run:        $(EXEPRJ)
	    @cd $(ANIBIN) && ./aniPRJ.exe

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.gmv 

help:		
	    @echo "make {exe|run|gmv-mesh1|gmv-mesh2|gmv-meta|gmv-remap|clean|help}"
	    @echo " "
	    @echo "   exe - create four executables in bin/"
	    @echo "   run - run the remap test"
	    @echo " "
	    @echo "   gmv-mesh1 - draw the first mesh"
	    @echo "   gmv-mesh2 - draw the second mesh"
	    @echo "   gmv-meta  - draw the metamesh"
	    @echo "   gmv-remap - draw the projected function"
	    @echo " "
	    @echo "   clean  - clean auto-built files"
	    @echo "   help   - print this help message"
	    @echo " "

info:	
	    @echo " " 
	    @echo "Executables are located in ../../../bin/"
	    @echo " "

gmv-mesh1:
	    @cd $(ANIBIN) && $(VIEWER) -a $(THISDIR)/try.attr -i mesh1.gmv

gmv-mesh2:
	    @cd $(ANIBIN) && $(VIEWER) -a $(THISDIR)/gmvM.attr -i mesh2.gmv

gmv-meta:
	    @cd $(ANIBIN) && $(VIEWER) -a $(THISDIR)/gmvM.attr -i metamesh.gmv

gmv-remap:
	    @cd $(ANIBIN) && $(VIEWER) -a $(THISDIR)/gmvU.attr -i remap.gmv

############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


$(EXEPRJ): $(OBJPRJ)
	@echo ''
	@echo 'Creating executable $(EXEPRJ)'
	@$(FLINKER) $(LDFLAGS) -o $(EXEPRJ) $(OBJPRJ) $(LIBS)
	@echo ''


