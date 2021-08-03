############################################################
ANIHOME = ../../..
include $(ANIHOME)/src/Rules.make

ANIBIN = $(ANIHOME)/bin

############################################################
FLAGS    := -W -Wall -O2 -I$(ANIINC)
CFLAGS   := $(CFLAGS) $(FLAGS)
CXXFLAGS := $(CXXFLAGS) $(FLAGS)
LDFLAGS  := $(LDFLAGS) -L$(ANILIB)
LDLIBS   := $(LDLIBS) -lm

ifdef CASROOT
-include $(ANILIB)/CGM/cgm.make
#export LD_LIBRARY_PATH=$(CASROOT)/Linux/lib
endif

############################################################
EXEAFT = $(ANIBIN)/FrontToTet_aft.exe
OBJAFT = main_aft.o

EXEPRM = $(ANIBIN)/ParamToTet_aft.exe
OBJPRM = main_prm.o

EXEMDF = $(ANIBIN)/FrtModifToTet_aft.exe
OBJMDF = main_mdf.o

EXESCG = $(ANIBIN)/ScgToTet_aft.exe
OBJSCG = main_scg.o

ifdef CASROOT
EXECAD = $(ANIBIN)/CadToTet_aft.exe
OBJCAD = main_cad.o
$(OBJCAD): CXXFLAGS:=$(CXXFLAGS) $(CGM_INCLUDES)
endif

############################################################
all: help

help:           
	@echo "make {exe|run-aft|run-prm|run-mdf|run-scg|run-cad|clean|gmv|help}"
	@echo " "
	@echo "   exe     - creates all executables"
	@echo "   run-aft - run executable"
	@echo "   run-prm - run executable"
	@echo "   run-mdf - run executable"
	@echo "   run-scg - run executable"
	@echo "   run-cad - run executable"
	@echo "   gmv     - visualize mesh using GMV"
	@echo "   test    - run batch tests"
	@echo "   clean   - clean the package"
	@echo "   help    - print this message"
	@echo " "


exe: $(EXEAFT) $(EXEPRM) $(EXEMDF) $(EXESCG) $(EXECAD)

run-aft: $(EXEAFT)
	 @cd $(ANIBIN) && ./FrontToTet_aft.exe ../data/rbox.frt 0.1
	 @echo "It is strongly recommended to postprocess the produced mesh"
	 @echo "bin/mesh.out by running bin/mbaFixAft.exe. To this end," 
	 @echo "cd ../PackageMBA; make exe; make run-fix; cd ../PackageAFT"

run-prm: $(EXEPRM)
	 @cd $(ANIBIN) && ./ParamToTet_aft.exe
	 @echo "It is strongly recommended to postprocess the produced mesh"
	 @echo "bin/mesh.out by running bin/mbaFixAft.exe. To this end," 
	 @echo "cd ../PackageMBA; make exe; make run-fix; cd ../PackageAFT"

run-mdf: $(EXEMDF)
#	 @cd $(ANIBIN) && ./FrtModifToTet_aft.exe ../data/rbox.frt            # for frt version
	 @cd $(ANIBIN) && ./FrtModifToTet_aft.exe ../data/v2_01_rescale.vtk   # for vtk version
	 @echo "It is strongly recommended to postprocess the produced mesh"
	 @echo "bin/mesh.out by running bin/mbaFixAft.exe. To this end," 
	 @echo "cd ../PackageMBA; make exe; make run-fix; cd ../PackageAFT"

run-scg: $(EXESCG)
	 @cd $(ANIBIN) && ./ScgToTet_aft.exe
	 @echo "It is strongly recommended to postprocess the produced mesh"
	 @echo "bin/mesh.out by running bin/mbaFixAft.exe. To this end," 
	 @echo "cd ../PackageMBA; make exe; make run-fix; cd ../PackageAFT"

ifdef CASROOT
run-cad: $(EXECAD)
	 @cd $(ANIBIN) && ./CadToTet_aft.exe ../data/rotate.brep
#	 @cd $(ANIBIN) && ./CadToTet_aft.exe ../data/box.brep
	 @echo "It is strongly recommended to postprocess the produced mesh"
	 @echo "bin/mesh.out by running bin/mbaFixAft.exe. To this end," 
	 @echo "cd ../PackageMBA; make exe; make run-fix; cd ../PackageAFT"

else
run-cad: $(EXECAD)
	 @echo '$$CASROOT is not set, OpenCASCADE support is disabled'


endif

clean:
	 rm -f *.o $(EXEAFT) $(EXEPRM) $(EXEMDF) $(EXESCG) $(EXECAD) \
	     $(ANIBIN)/mesh.??? $(ANIBIN)/frt.??? $(ANIBIN)/fail.??? $(ANIBIN)/refined.??? \
	     $(ANIBIN)/run_test.sh $(ANIBIN)/*.brep

gmv:        
	 @cd $(ANIBIN); $(VIEWER) -i mesh.gmv

test: $(EXEAFT)
	 cp -f run_test.sh $(ANIBIN)
	 #cd $(ANIBIN); for x in $(ANIDAT)/{cube,sculpture,rbox,knot,ball,test,gear1,gear6,gear10,cad1,cad2,cad3,hard,bumpy,dragon}.frt; do ./run_test.sh $$x; done
	 cd $(ANIBIN); for x in $(ANIDAT)/*.frt; do ./run_test.sh $$x; done


############################################################
$(EXEAFT):   $(OBJAFT)   
	$(LD) $(LDFLAGS) -o $(EXEAFT)  $(OBJAFT)  $(LIBAFT) $(LDLIBS)

$(EXEPRM):   $(OBJPRM)  
	$(LD) $(LDFLAGS) -o $(EXEPRM)  $(OBJPRM)  $(LIBAFT) $(LIBFRTPRM) $(LDLIBS)

$(EXEMDF):   $(OBJMDF)
	$(LD) $(LDFLAGS) -o $(EXEMDF)  $(OBJMDF)  $(LIBAFT) $(LIBFRTMDF) $(LDLIBS)
	
$(EXESCG):   $(OBJSCG) 
	$(LD) $(LDFLAGS) -o $(EXESCG)  $(OBJSCG)  $(LIBAFT) $(LIBFRTSCG) $(LIBFRTPRM) $(LDLIBS)

ifdef CASROOT
$(EXECAD):   $(OBJCAD) 
	$(LD) $(LDFLAGS) -o $(EXECAD)  $(OBJCAD)  $(LIBAFT) $(LIBFRTCAD) $(LIBFRTPRM) $(LDLIBS) $(CGM_LIBS_LINK) $(OCC_LIBS_LINK) -ldl

endif
