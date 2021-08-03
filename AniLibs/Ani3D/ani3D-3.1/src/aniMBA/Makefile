ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
OBJCORE	= mba_analytic.o mba_nodal.o mba_fixshape.o \
	  ani2.o makM.o check_mesh.o forlibmba.o \
	  makQ.o tangled.o list.o \
	  auxSP.o auxSR.o auxSF.o auxSE.o checks.o update.o \
	  tetrahedron.o triangle.o \
          F2E.o E2F.o swapR.o insrtP.o delete_point.o move_point.o collapse_edge.o split_tet.o \
          minim.o nlnfnc.o \
          lintrp3D.o dsort.o utils.o utils_add.o missing.o refine.o smoothing.o \
	  ZZ.o datacopy.o \
          error.o output.o isnan.o time.o time_new.o random.o

OBJIO	= loadM.o saveM.o 

OBJLIB	= $(OBJCORE) $(OBJIO)


############################################################
all:	    help

lib:        $(LIBMBA) info

clean:
	    @/bin/rm -f *.o *~ 

help:		
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib   - creates library libmba3D-$(version).a"
	    @echo "   clean - clean this directory"
	    @echo "   help  - print help message"
	    @echo " "

info:	
	    @echo " " 
	    @echo "Libriry version: " $(version)
	    @echo "Libraries are located in ../../lib/"
	    @echo " "


############################################################
.f.o:
	@echo '$(F77) -c $(FFLAGS) ' $*.f
	@$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	@$(CC)  $(CFLAGS) -c $*.c -o $*.o


$(LIBMBA): $(OBJLIB)
	@ar -r $(LIBMBA) $(OBJLIB)


############################################################
mba_analytic.o     : magic.fd  status.fd  output.fd
mba_nodal.o        : magic.fd  status.fd  output.fd
mba_fixshape.o     : magic.fd  status.fd  output.fd
metric3D.o         : magic.fd

ani2.o       	   : makS.fd  color.fd  magic.fd  output.fd  status.fd  operat.fd 

forlibfem.o        : metric.fd

saveM.o  	   : color.fd
makM.o   	   : makS.fd  color.fd  status.fd  magic.fd
check_mesh.o   	   : makS.fd  color.fd  status.fd 
makQ.o   	   : magic.fd cubature.fd          metric.fd

F2E.o    	   : makS.fd  color.fd  operat.fd  shutF.fd            status.fd
E2F.o    	   : makS.fd  color.fd  operat.fd  shutF.fd  magic.fd  status.fd
swapR.o  	   : makS.fd  color.fd  operat.fd  shutF.fd  magic.fd  status.fd
insrtP.o 	   : makS.fd  color.fd  operat.fd                      status.fd
deletP.o 	   : makS.fd  color.fd  operat.fd  magic.fd            status.fd
moveP.o  	   : makS.fd  color.fd  operat.fd  shutF.fd  magic.fd  status.fd
splitE.o 	   : makS.fd  color.fd  operat.fd

minim.o 	   : makS.fd  status.fd
nlnfnc.o	   : color.fd
auxSP.o 	   : makS.fd  color.fd  operat.fd
auxSR.o 	   : makS.fd  color.fd
auxSF.o 	   : magic.fd
auxSE.o 	   : color.fd
update.o 	   : makS.fd  color.fd

lintrp3D.o 	   : lintrp.fd
lintrp3D_test.o    : lintrp.fd

utils.o            : makS.fd
output.o           : output.fd

error.o  	   : output.fd
time.o   	   : time.c


