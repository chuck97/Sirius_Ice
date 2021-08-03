AFTHOME  := $(shell cd ../../..;pwd)
include $(AFTHOME)/src/cflags.make

############################################################
LIBSCG        = $(AFTLIB)/libfrtscg-$(version).a
OBJLIBSCG     = mesh.o      libfrtscg.o    libprim.o


############################################################
all: lib

lib: $(LIBSCG)

clean:
	rm -f *.o

purge: clean
	rm -f $(LIBSCG) 


############################################################
$(LIBSCG):    $(OBJLIBSCG)
	rm -f $(LIBSCG)
	$(AR) -r $(LIBSCG)  $(OBJLIBSCG)

