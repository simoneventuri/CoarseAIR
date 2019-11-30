#Fortran compiler given below
#FC=gfortran 

ifeq ($(FC),gfortran)
#FCFLAGS=-ffree-form -pipe -std=gnu -g -pg -O0 -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -Wsurprising -Wtabs -Wunderflow -Wintrinsic-shadow -ggdb -fdump-core -fno-unsafe-math-optimizations -finit-integer=n -finit-real=nan
FCFLAGS=-O0 -fno-unsafe-math-optimizations -Wmaybe-uninitialized
endif

DIR=$(PWD)
SRC=$(DIR)/src/
OBJ=$(DIR)/obj/
MOD=$(DIR)/mod/
EXEC=$(DIR)/exec/

#------------------
#Source files
src=utils_tecplot_binary.F90\
    ReCheckOriginalPESsAndDiatomic.F90

obj = $(addprefix $(OBJ), $(src:.F90=.o))

#------------------
#Includes 
INC=

#Modules
ifeq ($(FC),gfortran)
INC+=-J $(MOD) 
endif

#TecIO library
INC+=-I$(TECIO_INC)

LIBS=-L$(TECIO_LIB) -ltecio -lpthread

default: plot

$(OBJ)%.o: $(SRC)%.F90
	$(FC) $(INC) -c $(FCFLAGS) $< -o $@;


plot: $(obj)
	@echo 'Building plot program for N2O PES...'
	@$(FC) $(FCFLAGS) -o $(EXEC)plot_ $(obj) $(LIBS)
	@echo 'Done!'		

clean:
	@echo 'Cleaning...'
	@rm -f $(EXEC)box_
	@rm -f $(MOD)*.mod
	@rm -f $(OBJ)*.o
	@rm -f $(EXEC)*.out
	@echo 'Done!'
