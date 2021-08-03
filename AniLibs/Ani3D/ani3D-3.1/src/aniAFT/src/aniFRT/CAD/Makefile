AFTHOME  := $(shell cd ../../..;pwd)
include $(AFTHOME)/src/cflags.make

ifdef CASROOT
-include $(AFTLIB)/CGM/cgm.make
cgmwrap.o: CXXFLAGS := ${CXXFLAGS} ${CGM_INCLUDES} -DTEMPLATE_DEFS_INCLUDED

############################################################
LIBFRTCAD    = $(AFTLIB)/libfrtcad-$(version).a
OBJLIBFRTCAD = supportc.o    region3c.o    user3c.o \
               cgmwrap.o     libfrtcad.o

############################################################
all: lib

lib: $(LIBFRTCAD)

clean:
	rm -f *.o

purge: clean
	rm -f $(LIBFRTCAD)
	+cd $(AFTHOME)/../cgm && make -f Makefile-ani3D purge

############################################################
$(LIBFRTCAD): $(OBJLIBFRTCAD)
	rm -f $(LIBFRTCAD)
	$(AR) -r $(LIBFRTCAD)  $(OBJLIBFRTCAD)

cgmwrap.o: cgmwrap.cpp  cgm

cgm:
	+cd $(AFTHOME)/../cgm && make -f Makefile-ani3D lib

    
else
all: lib

lib:
	@echo '$$CASROOT is not set, OpenCASCADE support is disabled'

clean:

purge:

endif
