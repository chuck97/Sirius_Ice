############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make


############################################################
LIBS = -L$(ANILIB) $(LIBLU) $(LIBBLAS)


############################################################
EXEC = $(ANIBIN)/aniLU.exe


############################################################
OBJS = main.o aux_umf.o 


############################################################
all:    help

exe:    $(EXEC)

run:    $(EXEC)
	cd $(ANIBIN)  &&  ./aniLU.exe

help:
	@echo "make {exe|run|clean|help}"
	@echo " "
	@echo "   exe   - create one executable"
	@echo "   run   - run code with the matrix in bin/"
	@echo "   clean - clean this directory and bin"
	@echo " "

clean:
	rm -rf *.o  $(EXEC) $(ANIBIN)/input_main


############################################################
$(EXEC) : $(OBJS)
	  @$(F77) $(LFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	  @echo " "
	  @echo "Executables are located in ../../../bin"
	  @echo " "


############################################################
.c.o :
	$(CC)  $(CFLAGS) -c $*.c -o $*.o
.f.o :
	$(F77) $(FFLAGS) -c $*.f -o $*.o


############################################################
main.o	  	  : main.f 
aux_umf.o         : aux_umf.f

