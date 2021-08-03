ANIHOME = ./
#ANIHOME := $(shell /bin/pwd)
include $(ANIHOME)/src/Rules.make

########################################################################
TODAY   := $(shell date +%m-%d-%y)


########################################################################
DISTRO  = COPYRIGHT AUTHORS README DEMO ChangeLog ACKNOWLEDGEMENTS LICENSE VERSION \
	  bin lib include/* src/* cmake/* CMakeLists.txt Makefile \
 	  data/* doc/u*pdf 
 

########################################################################
all: 	help
weekly: clean snap


########################################################################
libs:
	@echo "======================================="
	@echo "        Making libraries               "
	@echo "======================================="
	@cd $(ANIHOME)/src/blas;    make lib
	@cd $(ANIHOME)/src/lapack;  make lib
	@cd $(ANIHOME)/src/cgm;     make -f Makefile-ani3D lib
	@cd $(ANIHOME)/src/umfpack; make lib
	@cd $(ANIHOME)/src/aniAFT;  make lib
	@cd $(ANIHOME)/src/aniFEM;  make lib
	@cd $(ANIHOME)/src/aniILU;  make lib
	@cd $(ANIHOME)/src/aniLU;   make lib
	@cd $(ANIHOME)/src/aniINB;  make lib
	@cd $(ANIHOME)/src/aniLMR;  make lib
	@cd $(ANIHOME)/src/aniMBA;  make lib
#	@cd $(ANIHOME)/src/aniMBAmpi; make lib
	@cd $(ANIHOME)/src/aniRCB;  make lib
	@cd $(ANIHOME)/src/aniPRJ;  make lib
	@cd $(ANIHOME)/src/aniVIEW; make lib
	@cd $(ANIHOME)/src/aniC2F;  make lib

clean:
	@dirlist="src/ani* src/lapack src/blas \
                src/Tutorials/Package* src/Tutorials/MultiPackage/[A,C,E,H,I,S]*"; \
        for dir in $$dirlist; do \
           cd $$dir; make clean; cd "$(CURDIR)"; \
        done; \
	cd src/cgm; make -f Makefile-ani3D clean; cd "$(CURDIR)";

distclean: 
	@cd $(ANIHOME)/src/cgm;      make -f Makefile-ani3D purge
	@cd $(ANIHOME)/src/aniAFT;   make purge
	@cd $(ANIHOME)/bin;   /bin/rm -fr *.exe

distro: $(DISTRO)
	tar -chvf  ani3D-3.1.tar --exclude='*MBAmpi' --exclude='*AMG' $(DISTRO) 
 
snap:   
	@echo "Preparing snapshot of the ani3d tree..."
	@tar -chvf dist/Ani3D-$(TODAY).tgz --gzip --exclude=*.pdf --exclude=*.tgz \
	     --exclude=*.bz2 --exclude='*.a' --exclude='*CVS*' --exclude='dist*'  \
	     --exclude=um*.tar.gz * > /dev/null

patch:	
	cd dist; \
	diff -ruNbB --exclude='CVS*'     --exclude='appl*'  --exclude='dist*'    \
	            --exclude='bin*'     --exclude='*.gmv'  --exclude='*.o'      \
	            --exclude='lib*'     --exclude='*.new'  --exclude='*.orig'   \
		    --exclude='*_dbg.f'  --exclude='*.pdf'  --exclude='doc*'     \
	            --exclude='*mpi*'    --exclude='data*'  --exclude='include*' \
		    --exclude='*CGMA*'   --exclude='testing*' \
	            ../../Ani3D-2.3 ani3D-2.3 > patch

help:
	@echo " "
	@echo "make {libs|clean|help|distro|weekly}"
	@echo " "



