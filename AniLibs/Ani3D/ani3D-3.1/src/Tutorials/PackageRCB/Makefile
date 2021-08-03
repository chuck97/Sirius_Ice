############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make
THISDIR := $(shell /bin/pwd)


############################################################
LIBS	= -L$(ANILIB) $(LIBVIEW) $(LIBRCB) $(LIBMBA) $(LIBBLAS) -lm 
	                

############################################################
EXERCB  = $(ANIBIN)/aniRCB.exe 


############################################################
RCBOBJS = main.o 


############################################################
all:	    help

exe:        $(EXERCB) info

clean:
	    @/bin/rm -f *.o *~ $(ANIBIN)/*.exe $(ANIBIN)/*.gmv

help:
	    @echo "make {exe|run|gmv-ini|gmv-rfn|gmv-crs|clean|help}"
	    @echo " "
	    @echo "   exe - compile the code"
	    @echo " "
	    @echo "   run - run the code"
	    @echo " "
	    @echo "   gmv-ini - draw initial mesh (bin/ini.gmv)"
	    @echo "   gmv-rfn - draw refined mesh (bin/rfn.gmv)"
	    @echo "   gmv-crs - draw coarsened mesh (bin/crs.gmv)"
	    @echo " "
	    @echo "   clean - clean the package"
	    @echo "   help  - print help message"
	    @echo " "

info:
	    @echo "Summary:"
	    @echo "Executables are located in bin/"
	    @echo " "

run:
	    @cd $(ANIBIN) && ./aniRCB.exe

gmv-ini:
	    @cd $(ANIBIN) && $(VIEWER) -a $(THISDIR)/gmv.attr -i ini.gmv

gmv-rfn:
	    @cd $(ANIBIN) && $(VIEWER) -a $(THISDIR)/gmv.attr -i rfn.gmv

gmv-crs:
	    @cd $(ANIBIN) && $(VIEWER) -a $(THISDIR)/gmv.attr -i crs.gmv


############################################################
.f.o:
	   @echo '$(F77) -c $(FFLAGS) ' $*.f
	   @$(F77) $(FFLAGS) -c $*.f -o $*.o 


############################################################
$(EXERCB): $(RCBOBJS) 
	   @echo ''
	   @echo 'Creating executable bin/aniRCB.exe'
	   @$(FLINKER) $(LDFLAGS) -o $(EXERCB) $(RCBOBJS) $(LIBS)


