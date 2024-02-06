#F90 = /sw/bin/g95
#F90 = /usr/bin/gfortran-4.2
#F77 = /sw/bin/fort77
BASE_PATH=/opt/intel/Compiler/11.1/072
MKL_INC=-I$(BASE_PATH)/mkl/include
F90COMP  = ifort
FFLAGS = -fpp  
FFTW_I_PATH = $(MKL_INC)/fftw 
#Executable extension:
EXT=.X
#Program name:
PRO=PFC
#executable target:
NAME=$(PRO)$(EXT)
 
#F90COMP  = gfortran
#FFLAGS = -x f95-cpp-input
#FFLAGS = 
########################################################################
MOD_OBJ =\
          modGlobalPFC.o modNoisePFC.o modCalcPFC.o  

F90_OBJS = pfcfft.o
##################################################################### Libraries

ALL_OBJS = $(MOD_OBJ) $(F90_OBJS) 

NDS_3D: $(ALL_OBJS)
	$(F90COMP) $(FFLAGS) $(ALL_OBJS) $(FFTW_I_PATH) -o PFC_BIN -mkl 

clean:  
	rm *.mod *.o PFC_BIN

#clean2:
#	rm out_* *.jpg
##################################################################### Rules
modGlobalPFC.o  : modGlobalPFC.f90
	          $(F90COMP) -c $(FFLAGS) $(FFTW_I_PATH) modGlobalPFC.f90 -mkl
modNoisePFC.o   : modGlobalPFC.o modNoisePFC.f90
	          $(F90COMP) -c $(FFLAGS)  $(FFTW_I_PATH) modNoisePFC.f90 -mkl
modCalcPFC.o    : modGlobalPFC.o modNoisePFC.o modCalcPFC.f90 
	          $(F90COMP) -c $(FFLAGS) $(FFTW_I_PATH) modCalcPFC.f90 -mkl
pfcfft.o   :  pfcfft.f90 modGlobalPFC.o modNoisePFC.o  modCalcPFC.o 
	          $(F90COMP) -c $(FFLAGS) $(FFTW_I_PATH) pfcfft.f90 -mkl
######################################################################cleanup 




