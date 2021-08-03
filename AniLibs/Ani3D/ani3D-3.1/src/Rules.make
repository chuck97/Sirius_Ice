############################################################
# Fortran and C compilers for SEQuential and PARallel codes.
############################################################
F77     = gfortran #compilers
CC      = gcc
CXX     = g++
VIEWER  = gmv      #visualization of GMV files

FLINKER = $(F77)   #linkers
CLINKER = $(CC)
LD      = $(CXX)  

LIBSYS  = -lgfortran        #for linking C and Fortran


############################################################
F77mpi      = /usr/lib64/openmpi/bin/mpif77
LINKERmpi   = $(F77mpi)
INCLUDE_MPI = -I /usr/include/openmpi-x86_64


############################################################
FFLAGS  = -O2
# FFLAGS  = -g -fbacktrace -ffpe-trap=zero,overflow,underflow 
CFLAGS  = -O2 
LDFLAGS = 


############################################################
version     = 3.0
UMFPACK_ver = 5.1
LAPACK_ver  = 3.2


ANIAFT	= $(ANIHOME)/src/aniAFT
ANIFEM 	= $(ANIHOME)/src/aniFEM
ANIILU  = $(ANIHOME)/src/aniILU 
ANIINB  = $(ANIHOME)/src/aniINB 
ANILMR  = $(ANIHOME)/src/aniLMR
ANIMBA	= $(ANIHOME)/src/aniMBA
ANIRCB	= $(ANIHOME)/src/aniRCB
ANIPRJ  = $(ANIHOME)/src/aniPRJ
ANIC2F  = $(ANIHOME)/src/aniC2F

ANIDOC 	= $(ANIHOME)/doc
ANILIB 	= $(ANIHOME)/lib
ANIBIN 	= $(ANIHOME)/bin
ANIDAT  = $(ANIHOME)/data
ANIINC  = $(ANIHOME)/include

ANIBLAS   = $(ANIHOME)/src/blas
ANILAPACK = $(ANIHOME)/src/lapack


############################################################
LIBAFT  = $(ANILIB)/libaft3D-$(version).a 
LIBFEM  = $(ANILIB)/libfem3D-$(version).a 
LIBILU  = $(ANILIB)/libilu-$(version).a 
LIBINB  = $(ANILIB)/libinb-$(version).a 
LIBLU   = $(ANILIB)/liblu-$(UMFPACK_ver).a 
LIBLMR  = $(ANILIB)/liblmr3D-$(version).a 
LIBMBA  = $(ANILIB)/libmba3D-$(version).a 
LIBRCB  = $(ANILIB)/librcb3D-$(version).a 
LIBPRJ  = $(ANILIB)/libprj3D-$(version).a 
LIBVIEW = $(ANILIB)/libview3D-$(version).a 
LIBC2F  = $(ANILIB)/libc2f3D-$(version).a

LIBFRTPRM = $(ANILIB)/libfrtprm-$(version).a
LIBFRTMDF = $(ANILIB)/libfrtmdf-$(version).a
LIBFRTSCG = $(ANILIB)/libfrtscg-$(version).a
LIBFRTCAD = $(ANILIB)/libfrtcad-$(version).a
# LIBBLAS   = -lblas
# LIBLAPACK = -llapack
LIBBLAS   = $(ANILIB)/libblas-$(LAPACK_ver).a
LIBLAPACK = $(ANILIB)/liblapack-$(LAPACK_ver).a $(ANILIB)/liblapack_ext-$(LAPACK_ver).a



############################################################
INCLUDE = -I$(ANIMBA) -I$(ANIFEM) -I$(ANIC2F)



