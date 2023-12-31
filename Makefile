# Default target
default : GITM

include Makefile.def
include Makefile.test

ABDIR   = srcSphereAB
MAINDIR = src
GLDIR   = srcGlow

PLANET=earth

src/ModSize.f90:
	cp src/ModSize.f90.orig src/ModSize.f90

install: src/ModSize.f90
	@(if [ ! -d srcData ]; then ln -s data/input srcData; fi)

# serial code uses NOMPI library
NOMPI:
	@echo "will make NOMPI"
	@echo ${NOMPIDIR}
	@cd ${NOMPIDIR}; make LIB

GITM:
	@cd ${SHAREDIR};           make LIB
	@cd $(ABDIR);              make -j1  LIB
	@cd $(EMPIRICALIEDIR);     make LIB
	@cd ${EMPIRICALUADIR};     make LIB
	@cd $(DATAREADINDICESDIR); make LIB
	@cd $(GLDIR);	           make -j1 LIB
	@cd $(MAINDIR);            make GITM

POST:
	@cd $(MAINDIR);  make POST

GITM2 = ${MYDIR}

LIB:
	cd $(ABDIR)     ; make                                         LIB
	cd $(GLDIR)     ; make LIBPREV=${GITM2}/${ABDIR}/libSphere.a   LIBADD
	cd $(MAINDIR)   ; make LIBPREV=${GITM2}/${GLDIR}/libUPTOGL.a   libGITM.a
	cd srcInterface ; make LIBPREV=${GITM2}/${MAINDIR}/libUA.a     LIB

serialrun:
	make GITM
	cd ${RUNDIR}; ${SERIAL} ./GITM.exe

clean:
	@cd $(ABDIR);    make clean
	@cd $(MAINDIR);  make clean
	@cd $(GLDIR);    make clean
	@cd srcInterface;make clean
	@(if [ -d share ]; then cd share; make clean; fi);
	@(if [ -d util ];  then cd util;  make clean; fi);

distclean: 
	./Config.pl -uninstall

allclean:
	@cd $(ABDIR);    make clean
	@cd $(MAINDIR);  make distclean
	@cd srcInterface;make distclean
	if [ -h srcData ]; then rm -f srcData; fi
#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/UA
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}; rm -f ./PostGITM.exe ; ln -s ${BINDIR}/PostProcess.exe ./PostGITM.exe
	cd ${RUNDIR}/UA; \
		mkdir restartOUT data DataIn; \
		ln -s restartOUT restartIN; \
		ln -s ${BINDIR}/pGITM .; \
		ln -s ${MYDIR}/srcData/* DataIn; rm -f DataIn/CVS
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/GITM.exe . ; \
		touch core ; chmod 444 core ; \
		ln -s UA/* .; \
	fi);

dist:
	make distclean
	tar cvzf gitm_`date "+%y%m%d"`.tgz Makefile* Config.pl get_info.pl \
	    share util src srcData srcDoc srcGlow srcIDL srcInterface \
	    srcPython srcMake srcSphereAB srcUser Copyright

