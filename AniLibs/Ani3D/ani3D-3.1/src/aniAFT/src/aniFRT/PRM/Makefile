AFTHOME  := $(shell cd ../../..;pwd)
include $(AFTHOME)/src/cflags.make

############################################################
LIBFRTPRM    = $(AFTLIB)/libfrtprm-$(version).a
OBJLIBFRTPRM = error3.o      libfrtprm.o   memory3.o \
	       refine3.o     region3.o     struct3.o \
	       support.o     tree32.o      tree3.o \
	       tria32.o      user3.o

############################################################
all: lib

lib: $(LIBFRTPRM)

clean:
	rm -f *.o

purge: clean
	rm -f $(LIBFRTPRM) 


############################################################
$(LIBFRTPRM):    $(OBJLIBFRTPRM)
	rm -f $(LIBFRTPRM)
	$(AR) -r $(LIBFRTPRM)  $(OBJLIBFRTPRM)

