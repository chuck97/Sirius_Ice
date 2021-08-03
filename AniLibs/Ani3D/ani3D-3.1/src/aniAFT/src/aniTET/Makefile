AFTHOME  := $(shell cd ../..;pwd)
include $(AFTHOME)/src/cflags.make

############################################################
LIBAFT     = $(AFTLIB)/libaft3D-$(version).a 
OBJLIBAFT  = aft2d.o    aft3d.o    helper.o \
             det.o      ugly.o     check.o \
             refine.o   libaft.o   delaunay.o \
	     cdt.o      qual.o

############################################################
all: lib

lib: $(LIBAFT)

clean:
	rm -f *.o

purge: clean
	rm -f $(LIBAFT)

############################################################
$(LIBAFT):   $(OBJLIBAFT) 
	rm -f $(LIBAFT)
	$(AR) -r $(LIBAFT)   $(OBJLIBAFT)

