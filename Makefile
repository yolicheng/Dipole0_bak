# This makefile works with the GNU make command, the one find on
# GNU/Linux systems and often called gmake on non-GNU systems, if you
# are using an old style make command, please see the file
# Makefile_oldstyle provided with the package.

# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
FCFLAGS = -O2

# the subfolders for modules, source files, object files, data, backup, doc-useful documents
SRC = ./src
MOD = ./mod
OBJ = ./obj

# absolute path for local libraries, specify for different computer
LIBDIR =  /home/yu/mylib  #/home/yu/mylib/

# flags to look for all modules (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I $(MOD)  #-I $(LIBDIR)/mod

# flags to put .mod files for compiled modules 
FMFLAGS = -J $(MOD)  

# libraries needed for linking,  -L declares the absolute path that contains the libraries to be used
# -llib states a library which is lib.a contained in this folder
#LDFLAGS = -li_need_this_lib

LDLAPACK = -L$(LIBDIR) -llapack -lblas   # general libraries, lapack and blas
LDEM = -L$(LIBDIR) -lemsys -ljpl -lemat  # used for polar in EM system

# List of executables to be built within the package
#PROGRAMS = prog1 prog2 prog3 prog4
#main = $(SRC)/vp1.f

#myobj = $(OBJ)/eigrg.o  $(OBJ)/dflrtz.o $(OBJ)/dflrtz_g3.o $(OBJ)/prt_eigval.o \
 #       $(OBJ)/prt_eigvec.o 

# myobj = $(OBJ)/$(wildcard *.o)
# myobj =  eigrg.o dflrtz.o dflrtz_g3.o prt_eigval.o prt_eigvec.o


#myobj = $(addprefix $(OBJ)/,  $(wildcard *.o) ) 
#   	  gr_rtbp.o  gr_cjrtbp.o  \ just for test
# gr_rk78.o   gr_lf.o gr_cjlf.o monomat.o  plob_fxd.o  plob_fxd.o   x2g.o invx2g.o  dflrtz_g3.o  dflrtz.o 

# update! 20160304 --- For the external subroutines the module calls, they should be complied before the module, or the module will be linked to empty external ones, without complain  and error.
# and use include to split the module into seperate files, avoiding tedious edit work. But it is strongly suggested that save all the related subroutines into one file. 

#dflrtz.o dflrtz_g3.o x2g.o invx2g.o gr_lf.o gr_cjlf.o-- all put into lfmod 

myobj = $(addprefix $(OBJ)/, gr_rk78.o plob.o inv.o  eigrg.o prt_eigval.o prt_eigvec.o \
 	 lfmod.o  pomod.o  monomat.o \
 	 gr_rtbp.o  gr_cjrtbp.o \
 	 gr_mpower.o  lf2cj.o \
 	 dtwofft.o dfour1.o prntft.o  plob_fxd.o fqmax.o fqext.o)

#myobj = $(addprefix $(OBJ)/, eigrg.o  dflrtz_g3.o prt_eigval.o prt_eigvec.o x2g.o invx2g.o \
#   	 lfmod.o  pomod.o  gr_lf.o gr_cjlf.o  dflrtz.o  gr_rtbp.o  gr_cjrtbp.o  \
# 	 monomat.o  gr_mpower.o mfdinst.o lf2cj.o  \
# 	 plob_fxd.o dtwofft.o dfour1.o prntft.o fqmax.o fqext.o)
# 	 
# poinc_n.o poinc_z0.o  --- to be modified
# join into pomod: pofam.o plob.o fctn.o dfcr.o deltx.o champ.o adams.o  sect.o poinc.o 
#myobj := $(addprefix $(OBJ)/,  $(patsubst %.c,%.o,$(wildcard *.c) ) )

#eigvp1
#vp_em: $(main) $(myobj)
#	$(FC) $(FCFLAGS) $(main) $(myobj) $(LDLAPACK) -o vp_em.exe 
	

mpo =  $(SRC)/main_po.f90
mpo1 =  $(SRC)/main_po1.f90
mpo2 =  $(SRC)/main_po2.f90

mpo21 =  $(SRC)/main_po21.f90
mpo22 =  $(SRC)/main_po22.f90
mpo23 =  $(SRC)/main_po23.f90
mpo24 =  $(SRC)/main_po24.f90
mpo25 =  $(SRC)/main_po25.f90
mpo26 =  $(SRC)/main_po26.f90
mpo27 =  $(SRC)/main_po27.f90
mpo28 =  $(SRC)/main_po28.f90


mpo3 =  $(SRC)/main_po3.f90

mpots =  $(SRC)/main_pots.f90

mmfd = $(SRC)/main_mfd.f90

meqmf = $(SRC)/main_eqmf.f90

#meq10 =  $(SRC)/main_eq10.f90
mpc10 =  $(SRC)/main_pc10.f90

mpc  =  $(SRC)/main_pc.f90

mfft =  $(SRC)/main_fft.f90
mfftst =  $(SRC)/main_fftst.f90

meq1 =  $(SRC)/main_eq1.f90

mfft1 =  $(SRC)/dxtwofft.f

mtest =  $(SRC)/test.f90

mtestphi =  $(SRC)/testphi.f90

testphi: $(mtestphi) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtestphi) $(myobj) $(LDLAPACK) -o testphi.exe 


test: $(mtest) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtest) $(myobj) $(LDLAPACK) -o test.exe 


pots: $(mpots) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpots) $(myobj) $(LDLAPACK) -o pots.exe 

po: $(mpo) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo) $(myobj) $(LDLAPACK) -o po.exe 

po1: $(mpo1) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo1) $(myobj) $(LDLAPACK) -o po1.exe 
	
po2: $(mpo2) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo2) $(myobj) $(LDLAPACK) -o po2.exe 
# to test the influence of the direction of the initial velocity, we test all the 8 cases, with \pm vz as one group, shown in one window
po21: $(mpo21) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo21) $(myobj) $(LDLAPACK) -o po21.exe	
po22: $(mpo22) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo22) $(myobj) $(LDLAPACK) -o po22.exe	
po23: $(mpo23) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo23) $(myobj) $(LDLAPACK) -o po23.exe	
po24: $(mpo24) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo24) $(myobj) $(LDLAPACK) -o po24.exe	
po25: $(mpo25) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo25) $(myobj) $(LDLAPACK) -o po25.exe	
po26: $(mpo26) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo26) $(myobj) $(LDLAPACK) -o po26.exe	
po27: $(mpo27) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo27) $(myobj) $(LDLAPACK) -o po27.exe	
po28: $(mpo28) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo28) $(myobj) $(LDLAPACK) -o po28.exe	
 

po3: $(mpo3) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo3) $(myobj) $(LDLAPACK) -o po3.exe 



mfd: $(mmfd) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mmfd) $(myobj) $(LDLAPACK) -o mfd.exe 

# eqmf  eqmf3 for bt=1 and bt=2 
eqmf: $(meqmf) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(meqmf) $(myobj) $(LDLAPACK) -o eqmf.exe 

pc: $(mpc) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpc) $(myobj) $(LDLAPACK) -o pc.exe 
 
pc10: $(mpc10) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpc10) $(myobj) $(LDLAPACK) -o pc10.exe 

fft: $(mfft) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfft) $(myobj) $(LDLAPACK) -o fft.exe 


fftst: $(mfftst) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfftst) $(myobj) $(LDLAPACK) -o fftst.exe 

eq1: $(meq1) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(meq1) $(myobj) $(LDLAPACK) -o eq1.exe 


fft1: $(mfft1) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfft1) $(myobj) $(LDLAPACK) -o fft1.exe 


# "make" builds all
#all: $(PROGRAMS)



   
# Another option is to put makefile in the SRC folder, but it then become complicated
# to transfer between different directories, for make and execute command
 
          
   
#apo: $(SRC)/main_apo.f90  $(myprogs) $(basicobj) 
#	$(FC) $(FCFLAGS) $(FMFLAGS) $(SRC)/main_apo.f90 $(myprogs) $(basicobj)  -o apo.exe 



#$intro
# ======================================================================
# Here comes the most interesting part: the rules for prog1, prog2,
# prog3 and prog4, modify to suit your needs
# ======================================================================

# In order to understand the next section, the process of building an
# executable has to be clear: this is typically done in two steps:

# 1.Compilation: every source file required for our program (tipically
# x.f90 or x.F90 in case of Fortran) is compiled into an object file
# (usually x.o)
# 2.Linking: the final executable file is built by "linking" together
# all the object files compiled in the previous step; in this step
# additional pre-compiled libraries of can be added, but this will not
# be treated here.

# These two steps are often performed through the same command and can
# be combined into a single operation, so it is easy to confuse them,
# but, in order to understand what comes further, one has to keep in
# mind that they are logically different operations.


#dependence
#prog2.o: prog2.incf
#prog3: aux.o
#prog4.o: mod.o
#prog4: mod.o

#$conclusion
# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) $(FMFLAGS) -o $@ -c $^ $(LDLAPACK)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
$(OBJ)/%.o: $(SRC)/%.f90
	$(FC) $(FCFLAGS) $(FMFLAGS) -o $@ -c $< 

$(OBJ)/%.o: $(SRC)/%.f
	$(FC) $(FCFLAGS) $(FMFLAGS) -o $@ -c $<


$(OBJ)/%.o: $(SRC)/%.F90
	$(FC) $(FCFLAGS) $(FMFLAGS) -o $@ -c $<
	

#$(OBJ)/%.o: $(SRC)/%.f90
#	$(FC) $(FCFLAGS) -c $<

#$(OBJ)/%.o: $(SRC)/%.f
#	$(FC) $(FCFLAGS) -c $<


#$(OBJ)/%.o: $(SRC)/%.F90
#	$(FC) $(FCFLAGS) -c $<
	

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f $(OBJ)/*.o  $(MOD)/*.mod  *.MOD *.o *.mod
first:
	mkdir src obj mod fig dat doc bak


veryclean: clean
	rm -f *~ $(PROGRAMS)


