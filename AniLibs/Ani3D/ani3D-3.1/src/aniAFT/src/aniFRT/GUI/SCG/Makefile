AFTHOME  := $(shell cd ../../../..;pwd)
AFTBIN := $(AFTHOME)/../../bin 
AFTLIB := $(AFTHOME)/../../lib 
include $(AFTHOME)/src/cflags.make

# correct LDLIBS/LDFLAGS/CFLAGS/CXXFLAGS if your freeglut package is installed in non-standard place

LD     := $(CXX)
LDLIBS := $(LDLIBS) -lglut -lGL
LDFLAGS:= $(LDFLAGS) -L$(AFTLIB) # -L/usr/local/lib
CFLAGS := $(CFLAGS) -I$(AFTHOME)/../../include # -I/usr/local/include
CXXFLAGS:=$(CXXFLAGS) -I$(AFTHOME)/../../include # -I/usr/local/include

############################################################
EXEGLUT       = $(AFTBIN)/gui_glut.exe
OBJGLUT       = gl2ps.o      gui_glut.o

LIBSCG        = $(AFTLIB)/libfrtscg-$(version).a
LIBAFT        = $(AFTLIB)/libaft3D-$(version).a
LIBFRTPRM     = $(AFTLIB)/libfrtprm-$(version).a


############################################################
all: exe

exe: $(EXEGLUT)

clean:
	rm -f *.o

purge: clean
	rm -f $(EXEGLUT)

############################################################
$(EXEGLUT):   $(OBJGLUT)  $(LIBSCG)  $(LIBFRTPRM)  $(LIBAFT)
	$(LD) $(LDFLAGS) -o $(EXEGLUT)  $(OBJGLUT)  $(LIBSCG)  $(LIBFRTPRM)  $(LIBAFT)  $(LDLIBS)

