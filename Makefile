#---------------------------------------------
#
#  Makefile for GNU make version 3.75 or later
#
#---------------------------------------------

BUILD = Tru64

COMPILETYPE = OPT
#COMPILETYPE = DEBUG
#COMPILETYPE = PROFILE


#
#  Alpha Tru64 Options
#

Tru64-F90COMPILER = f95
Tru64-DEBUGFLAGS = -g -C -warn argument_checking #-synchronous_exceptions -ladebug -DDEBUG
Tru64-OPTFLAGS = -g0 -fast -pipeline -transform_loops #-fast -O5 -speculate all -pipeline -transform_loops -non_shared
Tru64-PROFILEFLAGS = -p -g3 -fast #-pipeline -non_shared -transform_loops -speculate all -inline speed
Tru64-DEFINES = -DFFTW -DPNTSIZE=8 #-DMPI -DFFTW -DSSL2 -DCXML
Tru64-F90FLAGS = $(Tru64-$(COMPILETYPE)FLAGS) -cpp $(Tru64-DEFINES) -module ./obj/$(BUILD) #-convert little_endian 
Tru64-LIBS = -lfftw_mpi -lfftw -lmpi -lmpio -lm #-lVT -lpmpi -lfftw_mpi -lmpio -lmpi 
Tru64-LDFLAGS = -L/opt/fftw-2.1.3/lib #-L$$VAMPIR_ROOT/lib -om -WL,-om_dead_code -L/usr/local/lib -L/opt/fftw-2.1.3/lib
Tru64-INCLUDES = -I./obj/$(BUILD) #-I/usr/local/include
Tru64-MCC = -fast



#
#  MacOSX Options
#

MacOSX-F90COMPILER = f95
MacOSX-DEBUGFLAGS = -g #-C -warn argument_checking #-synchronous_exceptions -ladebug -DDEBUG
MacOSX-OPTFLAGS = #-g0 #-fast -O5 -speculate all -pipeline -transform_loops -non_shared
MacOSX-PROFILEFLAGS = #-p -g3 -fast #-pipeline -non_shared -transform_loops -speculate all -inline speed
MacOSX-DEFINES =-DFFTW -DPNTSIZE=4 #-DMPI -DFFTW -DSSL2 -DCXML
MacOSX-F90FLAGS = $(MacOSX-$(COMPILETYPE)FLAGS) $(MacOSX-DEFINES) -cpp -module ./obj/$(BUILD) #-convert little_endian 
MacOSX-LIBS = #-lfftw_mpi -lfftw -lmpi -lmpio -lm #-lVT -lpmpi -lfftw_mpi -lmpio -lmpi 
MacOSX-LDFLAGS = #-L/home/900/sdm900/lib #-L$$VAMPIR_ROOT/lib -om -WL,-om_dead_code -L/usr/local/lib -L/opt/fftw-2.1.3/lib
MacOSX-INCLUDES = -I./obj/$(BUILD) #-I/usr/local/include
MacOSX-MCC = 



#
#  Intel Linux Options
#

Linux-x86-F90COMPILER = ifc
Linux-x86-DEBUGFLAGS = -g -C
Linux-x86-OPTFLAGS = -O3 -tpp7 -xW
Linux-x86-LINKOPTFLAGS = 
Linux-x86-PROFILEFLAGS = 
Linux-x86-DEFINES = -DFFTW -DPNTSIZE=4
Linux-x86-F90FLAGS = $(Linux-x86-$(COMPILETYPE)FLAGS) -cpp $(Linux-x86-DEFINES) -module ./obj/$(BUILD) -I./obj/$(BUILD)
Linux-x86-LINKF90FLAGS = $(Linux-x86-LINKOPTFLAGS) $(Linux-x86-DEBUGFLAGS) $(Linux-x86-PROFILEFLAGS) $(Linux-x86-DEFINES)
Linux-x86-LIBS = -lfftw -lPEPCF90 -lIEPCF90
Linux-x86-LDFLAGS = -L/opt/fftw-2.1.4/lib -L/opt/intel/compiler70/ia32/lib
Linux-x86-INCLUDES = -I./obj/$(BUILD)
Linux-x86-MCC = 



#
#  Alpha Linux Options
#

Linux-alpha-F90COMPILER = fort
Linux-alpha-DEBUGFLAGS = $(OSF1-DEBUGFLAGS)
Linux-alpha-OPTFLAGS = $(OSF1-OPTFLAGS)
Linux-alpha-PROFILEFLAGS = $(OSF1-PROFILEFLAGS)
Linux-alpha-DEFINES = $(OSF1-DEFINES)
Linux-alpha-F90FLAGS = $(Linux-alpha-$(COMPILETYPE)FLAGS) $(Linux-alpha-DEFINES) -convert little_endian -fpp -module ./obj/$(BUILD)
Linux-alpha-LIBS = -lcpml -lfftw -lcxml
Linux-alpha-LDFLAGS = -L/usr/local/lib
Linux-alpha-INCLUDES = -I./obj/$(BUILD) -I/usr/local/include
Linux-alpha-MCC = -O4 -fast



#
#  Source Files
#

BINSRC = 			src/p.f90 \
				src/cps_setup.f90 \
				src/cps_propagate.f90 \
				src/cps_finish.f90


SRC =				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/interpreter.f90 \
				src/fft.f90 \
				src/conversion.f90 \
				src/misc.f90 \
				src/potential.f90 \
				src/setup.f90 \
				src/hamiltonian.f90 \
				src/propagate.f90



#
#  Test Source Files
#

TESTSRC =			src/precision.f90 \
				src/globals.f90 \
				src/fft.f90 \
				src/conversion.f90 \
				src/interpreter.f90 \
				src/misc.f90 \
				src/test.f90



#
#  Object Files
#

BINOBJECTS =  $(BINSRC:src/%.f90=obj/$(BUILD)/%.o)

OBJECTS =  $(SRC:src/%.f90=obj/$(BUILD)/%.o)



DIRS = 				bin \
				obj \
				bin/$(BUILD) \
				obj/$(BUILD)



#
#  Test Object Files
#

TESTOBJECTS = $(TESTSRC:src/%.f90=obj/$(BUILD)/%.o)



#
#  Some Extra Objects
#

EXTRAOBJS =



#
#  The actual binaries that you might want to build
#

p:				bin/$(BUILD)/p

cps:				bin/$(BUILD)/cps_setup \
				bin/$(BUILD)/cps_propagate \
				bin/$(BUILD)/cps_finish

test:				bin/$(BUILD)/test

binary:				bin/$(BUILD)/binary



#
#  The ugly stuff to build the binaries
#

bin/$(BUILD)/p: 		$(DIRS) \
				$(OBJECTS) \
				obj/$(BUILD)/p.o
	$($(BUILD)-F90COMPILER) $($(BUILD)-F90FLAGS) $($(BUILD)-INCLUDES) $(OBJECTS) obj/$(BUILD)/p.o $(EXTRAOBJS) -o $@ $($(BUILD)-LDFLAGS) $($(BUILD)-LIBS)



bin/$(BUILD)/cps_setup: 	$(DIRS) \
				$(OBJECTS) \
				obj/$(BUILD)/cps_setup.o
	$($(BUILD)-F90COMPILER) $($(BUILD)-F90FLAGS) $($(BUILD)-INCLUDES) $(OBJECTS) obj/$(BUILD)/cps_setup.o $(EXTRAOBJS) -o $@ $($(BUILD)-LDFLAGS) $($(BUILD)-LIBS)



bin/$(BUILD)/cps_propagate: 	$(DIRS) \
				$(OBJECTS) \
				obj/$(BUILD)/cps_propagate.o
	$($(BUILD)-F90COMPILER) $($(BUILD)-F90FLAGS) $($(BUILD)-INCLUDES) $(OBJECTS) obj/$(BUILD)/cps_propagate.o $(EXTRAOBJS) -o $@ $($(BUILD)-LDFLAGS) $($(BUILD)-LIBS)



bin/$(BUILD)/cps_finish: 	$(DIRS) \
				$(OBJECTS) \
				obj/$(BUILD)/cps_finish.o
	$($(BUILD)-F90COMPILER) $($(BUILD)-F90FLAGS) $($(BUILD)-INCLUDES) $(OBJECTS) obj/$(BUILD)/cps_finish.o $(EXTRAOBJS) -o $@ $($(BUILD)-LDFLAGS) $($(BUILD)-LIBS)



bin/$(BUILD)/test: 		$(DIRS) \
				$(TESTOBJECTS)
	$($(BUILD)-F90COMPILER) $($(BUILD)-F90FLAGS) $($(BUILD)-INCLUDES) $(TESTOBJECTS) -o $@ $($(BUILD)-LDFLAGS) $($(BUILD)-LIBS)



bin/$(BUILD)/binary: 	$(DIRS) \
				src/Makefile \
				src/binary.tm
	mcc $($(BUILD)-MCC) src/binary.tm -o $@



#
#  Object dependencies
#

obj/$(BUILD)/test.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/fft.f90 \
				src/conversion.f90 \
				src/interpreter.f90 \
				src/misc.f90 \
				src/test.f90 \
				src/Makefile



obj/$(BUILD)/p.o:		$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/misc.f90 \
				src/conversion.f90 \
				src/setup.f90 \
				src/hamiltonian.f90 \
				src/propagate.f90 \
				src/Makefile



obj/$(BUILD)/cps_setup.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/misc.f90 \
				src/conversion.f90 \
				src/setup.f90 \
				src/hamiltonian.f90 \
				src/propagate.f90 \
				src/Makefile



obj/$(BUILD)/cps_propagate.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/misc.f90 \
				src/conversion.f90 \
				src/setup.f90 \
				src/hamiltonian.f90 \
				src/propagate.f90 \
				src/Makefile



obj/$(BUILD)/cps_finish.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/misc.f90 \
				src/conversion.f90 \
				src/setup.f90 \
				src/hamiltonian.f90 \
				src/propagate.f90 \
				src/Makefile



obj/$(BUILD)/propagate.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/misc.f90 \
				src/conversion.f90 \
				src/hamiltonian.f90 \
				src/Makefile



obj/$(BUILD)/hamiltonian.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/fft.f90 \
				src/misc.f90 \
				src/potential.f90 \
				src/Makefile



obj/$(BUILD)/setup.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/fft.f90 \
				src/conversion.f90 \
				src/misc.f90 \
				src/Makefile



obj/$(BUILD)/potential.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90  \
				src/mpi.f90 \
				src/misc.f90 \
				src/interpreter.f90 \
				src/Makefile



obj/$(BUILD)/misc.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/interpreter.f90 \
				src/fft.f90 \
				src/conversion.f90 \
				src/Makefile



obj/$(BUILD)/conversion.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/Makefile



obj/$(BUILD)/fft.o:		$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/Makefile



obj/$(BUILD)/interpreter.o:	$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/mpi.f90 \
				src/Makefile



obj/$(BUILD)/mpi.o:		$(DIRS) \
				src/precision.f90 \
				src/globals.f90 \
				src/Makefile



obj/$(BUILD)/globals.o:	$(DIRS) \
				src/precision.f90 \
				src/Makefile



obj/$(BUILD)/precision.o:	$(DIRS) \
				src/Makefile



#
#  Build the actual object files
#

obj/$(BUILD)/%.o : 		src/%.f90
	$($(BUILD)-F90COMPILER) -c $($(BUILD)-F90FLAGS) $< -o $@



#
#  To create the appropriate directories
#

$(DIRS) :
	if ( test ! -d $@ ) then ( mkdir $@ ) fi



#
#  Clean the distribution
#

clean:
	rm -f obj/$(BUILD)/*
	rm -f bin/$(BUILD)/p
	rm -f bin/$(BUILD)/cps_*

clean-test:
	rm -f  bin/$(BUILD)/test
