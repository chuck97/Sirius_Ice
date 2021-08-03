ANIHOME = ../..
include $(ANIHOME)/src/Rules.make

############################################################
OBJLIB	= Nodal2MetricVAR.o   Nodal2MetricZZ.o \
	  Func2MetricMAX.o    Func2GradMetricMAX.o \
	  EdgeEst2MetricMAX.o EdgeEst2MetricLS.o EdgeEst2GradMetricMAX.o \
	  CellEst2MetricZZ.o  Lp_norm.o utils.o error.o


############################################################
all:	    help

lib:        $(LIBLMR) info

clean:
	    @/bin/rm -f *.o *~ 

help:		
	    @echo "make {lib|clean|help}"
	    @echo " "
	    @echo "   lib   - creates library liblmr3D-$(version).a"
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


$(LIBLMR): $(OBJLIB)
	@ar -r $(LIBLMR) $(OBJLIB)


############################################################
metric_nodal.o     : metric.fd

