AM_CPPFLAGS += $(HAVE_OCC_DEF) \
               -I$(top_srcdir)/util \
               -I$(top_builddir)/util \
               -I$(top_srcdir)/geom \
               $(OCC_INC_FLAG)

# The name of the library to build
noinst_LTLIBRARIES = libcubit_OCC.la

libcubit_OCC_la_LIBADD = -L$(OCC_LIB_DIR) $(OCC_LIBS)

# The directory where headers will be installed
libcubit_OCC_la_includedir = $(includedir)

# The non-template sources
libcubit_OCC_la_SOURCES = \
    OCCShapeAttributeSet.cpp \
    OCCAttribSet.cpp \
    OCCBody.cpp \
    OCCCoEdge.cpp \
    OCCCurve.cpp \
    OCCLoop.cpp \
    OCCLump.cpp \
    OCCModifyEngine.cpp \
    OCCPoint.cpp \
    OCCQueryEngine.cpp \
    OCCShell.cpp \
    OCCSurface.cpp \
    OCCDrawTool.cpp

# Headers to be installed.  If any file in this list should
# not be installed, move it to the _SOURCES list above.
libcubit_OCC_la_include_HEADERS = \
    OCCShapeAttributeSet.hpp \
    OCCAttribSet.hpp \
    OCCBody.hpp \
    OCCCoEdge.hpp \
    OCCCurve.hpp \
    OCCGeometryCreator.hpp \
    OCCLoop.hpp \
    OCCLump.hpp \
    OCCModifyEngine.hpp \
    OCCPoint.hpp \
    OCCQueryEngine.hpp \
    OCCShell.hpp \
    OCCSurface.hpp \
    OCCDrawTool.hpp

