FC     	= ifort
LC     	= $(FC)
EXE    	= ngpdf
FITSDIR = #/home/wilmar/usr/local/lib/lib
LIBFITS = #cfitsio
INCDIR	= /opt/intel/intel/composer_xe_2015.0.090/mkl/include
IDIR	= /home/cardonac/projects/NG-PDF/ranlib/lib
LIBDIR	= /opt/intel/intel/composer_xe_2015.0.090/mkl/lib/intel64	
LDIR	= /opt/intel/intel/composer_xe_2015.0.090/compiler/lib/intel64
F_FL   	= -g -check all -mkl $(INCDIR)/mkl_dfti.f90 -openmp -I$(IDIR) #-DGFORTRAN -fno-second-underscore -fopenmp -fPIC -g
LIB_FL 	= -L$(LIBDIR) -mkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -L$(LDIR) -liomp5 

#####################
OBJ   =  ngpdf.o arrays.o fiducial.o functions.o 

def:	$(OBJ) $(OBJNR) $(OBJODE)
	$(LC) $(F_FL) $(OBJ) $(OBJNR) $(OBJODE) -o $(EXE)  $(LIB_FL)

%.o:	%.f90
	$(FC) $(F_FL) -c $<

%.o:	%.F90
	$(FC) $(F_FL) -c $<

%.o:	%.f
	$(FC) $(F_FL) -c $<

clean :
	rm -f *.o *.mod *.ini *~  fort.* slurm*.out $(EXE)

### put dependencies here ###

ngpdf.o :	arrays.o fiducial.o functions.o

