AFTHOME  := $(shell cd ../../..;pwd)
include $(AFTHOME)/src/cflags.make

############################################################
LIBFRTMDF    = $(AFTLIB)/libfrtmdf-$(version).a
OBJLIBFRTMDF = support.o     surfrf.o      3Dto2D1.o \
               polyreg.o     libfrtmdf.o

############################################################
all: lib

lib: $(LIBFRTMDF)

clean:
	rm -f *.o

purge: clean
	rm -f $(LIBFRTMDF)


############################################################
$(LIBFRTMDF):    $(OBJLIBFRTMDF)
	rm -f $(LIBFRTMDF)
	$(AR) -r $(LIBFRTMDF)  $(OBJLIBFRTMDF)

