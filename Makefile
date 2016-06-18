# Stuff for you; GNU
#FC = f95
#CC = gcc
#FC_PAR = /public/apps/openmpi/1.6.4/gnu.4.7.2/bin/mpif90 -DGNU
#CC_PAR = /public/apps/openmpi/1.6.4/gnu.4.7.2/bin/mpicc
# Options
#FFLAGS = -O2 -Wall -Wconversion -fbounds-check -Wuninitialized -fbacktrace
#FCFLAGS = $(FFLAGS)
#CFLAGS = -O2 -fbounds-check -Wuninitialized
#PLAT = LINUX
#ARCH = ./Obj
# Library directories
#VISITDIR = /home/bakerb3/libs
#MUMPSROOT = /home/bakerb3/libs
#METISDIR = /home/bakerb3/libs
#SCALAPDIR = /home/bakerb3/libs
#ITERDIR = /home/bakerb3/libs
# Stuff for you; INTEL
FC = /public/apps/intel/12.1.0/composer_xe_2011_sp1.6.233/bin/intel64/ifort
CC = /public/apps/intel/12.1.0/composer_xe_2011_sp1.6.233/bin/intel64/icc
FC_PAR = /public/apps/openmpi/1.6.4/intel.12.1.0/bin/mpif90 -DINTEL
CC_PAR = /public/apps/openmpi/1.6.4/intel.12.1.0/bin/mpicc
# Options
FFLAGS = -O2 -fpe0 -check uninit -check bounds -check pointers -align -ftrapuv -warn unused -assume byterecl -traceback
#FFLAGS = -O2 -check uninit -check pointers -check noarg_temp_created -traceback -fpe0 -align -cxxlib -ftrapuv -warn unused -assume byterecl
#FFLAGS = -O2 -traceback -fpe0 -align -cxxlib -ftrpuv -warn unused -assume byterecl
FCFLAGS = $(FFLAGS)
CFLAGS = -O2
#PLAT = LINUX
ARCH = ./Obj
# Library directories
VISITDIR = /home/bakerb3/intel/libs
MUMPSROOT = /home/bakerb3/intel/libs
METISDIR = /home/bakerb3/intel/libs
#SCALAPDIR = /home/ben/intel/scalapack-2.0.2
# Scalapack is tricky with MKL
INTELLIB = /public/apps/intel/12.1.0/composer_xe_2011_sp1.6.233/mkl/lib/intel64/lib/
LIBINTEL = -L$(INTELLIB)
LPBLAS = $(LIBINTEL) -lmkl_scalapack_lp64 -Wl,--start-group  \
	-lmkl_intel_lp64 -lmkl_sequential \
	-lmkl_core -lmkl_blacs_openmpi_lp64 -Wl,--end-group

#LIBLAPACK = /home/bakerb3/libs/liblapack.a /home/bakerb3/libs/libblas.a
LOTHERS = -lpthread -lm

# Orderings
#INCMUMPS = -I/home/bakerb3/include
INCMUMPS = -I$(MUMPSROOT)/
#INCVISIT = -I/home/bakerb3/include
INCVISIT = -I$(VISITDIR)
LMETIS = -L$(METISDIR) -lparmetis -lmetis
LPORD = -L$(MUMPSROOT)/PORD/lib -lpord
LORDERINGS = $(LSCOTCH) $(LMETIS) $(LPORD)
# Scalapack
#LIBSCALAP = -L$(SCALAPDIR) -lscalapack
#LPBLAS = $(LIBSCALAP) $(LIBLAPACK)
# MUMPS
LIBMUMPS = -L$(MUMPSROOT)/lib -lcmumps$(PLAT) -lmumps_common$(PLAT)
# VisIt visualization for paraview
LIBVISIT = -L$(VISITDIR) -lvisit
#INCVISIT = -I$(VISITDIR)
LIBALL = $(LIBMUMPS) $(LORDERINGS) $(LSCALAP) $(LPBLAS) $(LIBVISIT) $(LIBOTHERS)
INCALL = -I./ $(INCMUMPS) $(INCSCOTCH) $(INCVISIT)

XWIGGLE = xwiggle25
XWAVE25 = xbielak25
XMOVIE25 = xmovie
#XLBFGS = xlbfgs_hsl
#XDISP25 = xdisp25
XEXACT = xexact25
XGN = xgn25
XSTFREC = xsrcrec25
XJOINT = xjoint
EXECS = $(XWAVE25) $(XWIGGLE) $(XMOVIE25) $(XLBFGS) $(XDISP25) $(XEXACT) $(XGN) \
	$(XSTFREC) $(XJOINT)

BielakDriver = $(ARCH)/bielak25.o $(ARCH)/mesh_utils.o
FEMsrc = $(ARCH)/asmble.o $(ARCH)/csf2d.o $(ARCH)/bielak_utils.o \
	 $(ARCH)/dshl.o $(ARCH)/gaussq.o $(ARCH)/estiff25.o $(ARCH)/cjac.o \
	 $(ARCH)/locrec2d.o $(ARCH)/srclist.o $(ARCH)/fill1d.o
GraphUtils = $(ARCH)/graph25_driver.o $(ARCH)/graph25.o
MPIUtils = $(ARCH)/mkMPIgroups.o $(ARCH)/bcastfwd.o $(ARCH)/bcastinv.o
CommonUtils = $(ARCH)/sutils.o $(ARCH)/sort.o $(ARCH)/dft.o $(ARCH)/distaz.o \
	      $(ARCH)/lisdir.o $(ARCH)/wrapph.o $(ARCH)/window.o $(ARCH)/iir.o \
	      $(ARCH)/read_ini.o
SurfWaves = $(ARCH)/grns_driver.o $(ARCH)/grns_swr.o $(ARCH)/mkerthmods.o $(ARCH)/dfftl.o \
	    $(ARCH)/earthsr.o $(ARCH)/earthsubs.o \
	    $(ARCH)/srgramf.o $(ARCH)/mkhomogsrf.o
HaskSrc = $(ARCH)/haskgrn.o $(ARCH)/haskgrnsh.o $(ARCH)/haskinf.o \
	  $(ARCH)/haskattn.o $(ARCH)/genex.o
IOUtils = $(ARCH)/packdata.o $(ARCH)/vtk25.o $(ARCH)/dataio.o \
	  $(ARCH)/srcsub.o $(ARCH)/seismio.o $(ARCH)/invvtk25.o
SrcRec = $(ARCH)/srcrec.o $(ARCH)/mesh_utils.o $(ARCH)/srcupd.o \
	 $(ARCH)/recupd.o $(ARCH)/objfn.o $(ARCH)/jinv_utils.o \
	 $(ARCH)/adjgrad.o $(ARCH)/gnutils.o $(ARCH)/hybrd1.o
WigSrc = $(ARCH)/wiggle.o
LbfgsHsl = $(ARCH)/lbfgs_hsl.o $(ARCH)/lbfgs.o $(ARCH)/mesh_utils.o
GaussNewt = $(ARCH)/gn25.o $(ARCH)/mcsrch.o $(ARCH)/mesh_utils.o
JointDriver = $(ARCH)/sbgninv.o $(ARCH)/mcsrch.o $(ARCH)/mesh_utils.o \
	      $(ARCH)/fgh_window.o
InvUtils = $(ARCH)/djac.o $(ARCH)/dstiff25.o $(ARCH)/vforce25.o \
	   $(ARCH)/objfn.o $(ARCH)/adjgrad.o \
	   $(ARCH)/hessloc.o $(ARCH)/gengmask.o $(ARCH)/funcgradh.o \
	   $(ARCH)/precongrad.o $(ARCH)/pchadj.o \
	   $(ARCH)/gnmfree.o $(ARCH)/srcupd.o $(ARCH)/recupd.o \
	   $(ARCH)/qpmfree.o $(ARCH)/gnutils.o $(ARCH)/hybrd1.o \
	   $(ARCH)/jinv_utils.o
# mesh_utils.f bielak_utils.f
#FEMdriver = test.f90 mkmesh.f bielak.f
MovieDriver = $(ARCH)/vtkmovie.o $(ARCH)/mesh_utils.o
#MPIUtils = mkMPIgroups.f bcastfwd.f
#FEMsrc = graph_loc.f gensrc.f csf2d.f dshl.f gaussq.f cjakcm.f estiff.f \
	 haskgrn.f haskinf.f locate.f asmble.f90 asmble.f
#CommonUtils = dft.f90 sutils.f sort.f vtk2d.f dataio.f surfio.f \
	packdata.f
Disp25 = $(ARCH)/disp25.o $(ARCH)/graph25.o $(ARCH)/sort.o  \
	 $(ARCH)/csf2d.o    $(ARCH)/dshl.o    $(ARCH)/gaussq.o \
	 $(ARCH)cjac.o $(ARCH)/estiff25.o
Exact25 = $(ARCH)/exact.o $(ARCH)/mesh_utils.o $(ARCH)/dshl.o \
	  $(ARCH)/csf2d.o $(ARCH)/gaussq.o $(ARCH)/srclist.o   \
	  $(ARCH)/bielak_utils.o $(ARCH)/locrec2d.o $(ARCH)/fill1d.o
Csrc = $(ARCH)/visit_wrapper.o

OBJ1 = $(BielakDriver) \
	$(CommonUtils) \
	$(GraphUtils) \
	$(MPIUtils) \
	$(HaskSrc) \
	$(IOUtils) \
	$(FEMsrc) \
	$(SurfWaves) \
	$(Csrc)

OBJ2 = $(WigSrc) \
	$(CommonUtils) \
	$(IOUtils) \
	$(Csrc)

OBJ3 = $(MovieDriver) \
	$(CommonUtils) \
	$(IOUtils) \
	$(Csrc)

OBJ4 = $(LbfgsHsl) \
	$(CommonUtils) \
	$(GraphUtils) \
	$(MPIUtils) \
	$(FEMsrc) \
	$(HaskSrc) \
	$(InvUtils) \
	$(IOUtils) \
	$(SurfWaves) \
	$(Csrc)

OBJ5 = $(Disp25)

OBJ6 = $(Exact25) \
	$(HaskSrc) \
	$(SurfWaves) \
	$(CommonUtils) \
	$(GraphUtils) \
	$(IOUtils) \
	$(MPIUtils) \
	$(Csrc)

OBJ7 =  $(GaussNewt) \
	$(CommonUtils) \
	$(GraphUtils) \
	$(MPIUtils) \
	$(FEMsrc) \
	$(HaskSrc) \
	$(InvUtils) \
	$(IOUtils) \
	$(SurfWaves) \
	$(Csrc)

OBJ8 =  $(SrcRec) \
	$(CommonUtils) \
	$(GraphUtils) \
	$(MPIUtils) \
	$(HaskSrc) \
	$(IOUtils) \
	$(FEMsrc) \
	$(SurfWaves) \
	$(Csrc)

OBJ9  = $(JointDriver) \
	$(CommonUtils) \
	$(GraphUtils) \
	$(MPIUtils) \
	$(FEMsrc) \
	$(HaskSrc) \
	$(InvUtils) \
	$(IOUtils) \
	$(SurfWaves) \
	$(Csrc)

all: $(EXECS)

$(XWAVE25): $(OBJ1)
	$(FC_PAR) $(FCFLAGS) -o $(XWAVE25) $(OBJ1) $(LIBALL) -I$(INCMUMPS)

$(XWIGGLE): $(OBJ2)
	$(FC_PAR) $(FCFLAGS) -o $(XWIGGLE) $(OBJ2) $(LIBVISIT)

$(XMOVIE25): $(OBJ3)
	$(FC_PAR) $(FCFLAGS) -o $(XMOVIE25) $(OBJ3) $(LIBVISIT)

$(XLBFGS): $(OBJ4)
	$(FC_PAR) $(FCFLAGS) -o $(XLBFGS) $(OBJ4) $(LIBALL)

#$(XDISP25): $(OBJ5)
#	$(FC_PAR) $(FCFLAGS) -o $(XDISP25) $(OBJ5) $(LMETIS) $(LIBLAPACK)

$(XEXACT): $(OBJ6)
	$(FC_PAR) $(FCFLAGS) -o $(XEXACT) $(OBJ6) $(LIBALL)

$(XGN): $(OBJ7)
	$(FC_PAR) $(FCFLAGS) -o $(XGN) $(OBJ7) $(LIBALL)

$(XSTFREC): $(OBJ8)
	$(FC_PAR) $(FCFLAGS) -o $(XSTFREC) $(OBJ8) $(LIBALL)

$(XJOINT): $(OBJ9)
	$(FC_PAR) $(FCFLAGS) -o $(XJOINT) $(OBJ9) $(LIBALL)

$(ARCH)/adjgrad.o: adjgrad.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/adjgrad.o adjgrad.f90

$(ARCH)/asmble.o: asmble.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/asmble.o asmble.f90

$(ARCH)/bielak25.o: bielak25.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/bielak25.o  bielak25.f90

$(ARCH)/dft.o: dft.f90
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/dft.o dft.f90

#$(ARCH)/disp25.o: disp25.f90 $(INCMUMPS)/cmumps_struc.h ./fwd_struc.h
#	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/disp25.f90

$(ARCH)/exact.o: exact.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/exact.o exact.f90

$(ARCH)/fill1d.o: fill1d.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/fill1d.o fill1d.f90

$(ARCH)/fgh_window.o: fgh_window.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/fgh_window.o fgh_window.f90

$(ARCH)/funcgradh.o: funcgradh.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/funcgradh.o funcgradh.f90

$(ARCH)/gengmask.o: gengmask.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/gengmask.o gengmask.f90

$(ARCH)/gn25.o: gn25.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/gn25.o gn25.f90

$(ARCH)/gnmfree.o: gnmfree.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/gnmfree.o gnmfree.f90

$(ARCH)/sbgninv.o: sbgninv.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/sbgninv.o sbgninv.f90

$(ARCH)/gnutils.o: gnutils.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/gnutils.o gnutils.f90

$(ARCH)/graph25_driver.o: graph25_driver.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/graph25_driver.o graph25_driver.f90

$(ARCH)/grns_driver.o: grns_driver.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/grns_driver.o grns_driver.f90

$(ARCH)/grns_swr.o: grns_swr.f90
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/grns_swr.o grns_swr.f90

$(ARCH)/hessloc.o: hessloc.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/hessloc.o hessloc.f90

$(ARCH)/hybrd1.o: hybrd1.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/hybrd1.o hybrd1.f

$(ARCH)/iir.o: iir.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/iir.o iir.f

$(ARCH)/jinv_utils.o: jinv_utils.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/jinv_utils.o jinv_utils.f90

$(ARCH)/qpmfree.o: qpmfree.f90
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/qpmfree.o qpmfree.f90

$(ARCH)/lbfgs_hsl.o: lbfgs_hsl.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/lbfgs_hsl.o lbfgs_hsl.f90

#$(ARCH)/migrate.o: migrate.f90 $(INCMUMPS)/cmumps_struc.h ./fwd_struc.h
#	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/migrate.o migrate.f90

$(ARCH)/mkerthmods.o: mkerthmods.f90
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/mkerthmods.o  mkerthmods.f90

$(ARCH)/pchadj.o: pchadj.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/pchadj.o pchadj.f90

$(ARCH)/read_ini.o: read_ini.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/read_ini.o read_ini.f90

$(ARCH)/srcrec.o: srcrec.f90 fwd_struc.h
	$(FC_PAR) $(FCFLAGS) $(INCMUMPS) -c -o $(ARCH)/srcrec.o srcrec.f90

$(ARCH)/vtkmovie.o: vtkmovie.f90
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/vtkmovie.o vtkmovie.f90

$(ARCH)/wiggle.o: wiggle.f90
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/wiggle.o wiggle.f90

$(ARCH)/bcastfwd.o: bcastfwd.f fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/bcastfwd.o bcastfwd.f

$(ARCH)/bcastinv.o: bcastinv.f fwd_struc.h
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/bcastinv.o bcastinv.f

$(ARCH)/bielak_utils.o: bielak_utils.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/bielak_utils.o bielak_utils.f

$(ARCH)/cjac.o: cjac.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/cjac.o cjac.f

$(ARCH)/csf2d.o: csf2d.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/csf2d.o csf2d.f

$(ARCH)/dataio.o: dataio.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/dataio.o dataio.f

$(ARCH)/dfftl.o: dfftl.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/dfftl.o dfftl.f

#$(ARCH)/dgamma.o: dgamma.f
#	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/dgamma.o dgamma.f

$(ARCH)/distaz.o: distaz.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/distaz.o distaz.f

$(ARCH)/djac.o: djac.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/djac.o djac.f

$(ARCH)/dshl.o: dshl.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/dshl.o dshl.f

$(ARCH)/dstiff25.o: dstiff25.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/dstiff25.o dstiff25.f

$(ARCH)/earthsr.o: earthsr.f ./commons.inc ./sizes.inc
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/earthsr.o earthsr.f

$(ARCH)/earthsubs.o: earthsubs.f ./commons.inc ./sizes.inc
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/earthsubs.o earthsubs.f
 
$(ARCH)/estiff25.o: estiff25.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/estiff25.o estiff25.f

$(ARCH)/genex.o: genex.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/genex.o genex.f

$(ARCH)/gaussq.o: gaussq.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/gaussq.o gaussq.f

$(ARCH)/gengrns.o: gengrns.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/gengrns.o gengrns.f

$(ARCH)/graph25.o: graph25.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/graph25.o graph25.f

$(ARCH)/haskattn.o: haskattn.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/haskattn.o haskattn.f

$(ARCH)/haskgrn.o: haskgrn.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/haskgrn.o haskgrn.f

$(ARCH)/haskgrnsh.o: haskgrnsh.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/haskgrnsh.o  haskgrnsh.f

$(ARCH)/haskinf.o: haskinf.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/haskinf.o haskinf.f

$(ARCH)/lbfgs.o: lbfgs.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/lbfgs.o lbfgs.f

$(ARCH)/invvtk25.o: invvtk25.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/invvtk25.o invvtk25.f

$(ARCH)/lisdir.o: lisdir.F
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/lisdir.o lisdir.F

$(ARCH)/locrec2d.o: locrec2d.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/locrec2d.o locrec2d.f

$(ARCH)/mcsrch.o: mcsrch.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/mcsrch.o mcsrch.f

$(ARCH)/mesh_utils.o: mesh_utils.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/mesh_utils.o mesh_utils.f

$(ARCH)/mkhomogsrf.o: mkhomogsrf.f ./commons.inc ./sizes.inc
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/mkhomogsrf.o mkhomogsrf.f

$(ARCH)/mkMPIgroups.o: mkMPIgroups.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/mkMPIgroups.o mkMPIgroups.f

$(ARCH)/objfn.o: objfn.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/objfn.o objfn.f

$(ARCH)/packdata.o: packdata.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/packdata.o packdata.f

$(ARCH)/precongrad.o: precongrad.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/precongrad.o precongrad.f

$(ARCH)/quadeig.o: quadeig.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/quadeig.o quadeig.f

$(ARCH)/recupd.o: recupd.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/recupd.o recupd.f

$(ARCH)/seismio.o: seismio.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/seismio.o seismio.f

$(ARCH)/sort.o: sort.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/sort.o sort.f

$(ARCH)/srclist.o: srclist.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/srclist.o srclist.f

$(ARCH)/srcsub.o: srcsub.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/srcsub.o srcsub.f

$(ARCH)/srcupd.o: srcupd.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/srcupd.o srcupd.f

$(ARCH)/srgramf.o: srgramf.f ./commons.inc ./sizes.inc
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/srgramf.o srgramf.f

$(ARCH)/sutils.o: sutils.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/sutils.o sutils.f

$(ARCH)/vforce25.o: vforce25.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/vforce25.o vforce25.f

$(ARCH)/vtk25.o: vtk25.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/vtk25.o vtk25.f

$(ARCH)/visit_wrapper.o: visit_wrapper.c
	$(CC_PAR) $(CFLAGS) $(INCVISIT) -c -o $(ARCH)/visit_wrapper.o visit_wrapper.c

$(ARCH)/window.o: window.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/window.o window.f

$(ARCH)/wrapph.o: wrapph.f
	$(FC_PAR) $(FCFLAGS) -c -o $(ARCH)/wrapph.o wrapph.f

clean:
	@$(RM) $(ARCH)/*.o $(EXECS)
