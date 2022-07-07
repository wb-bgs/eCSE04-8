main=gd2gc_cmdln

# decide which compiler to use the library depending on the architecture
FC = ftn
FFLAGS = -O3

#--------------------------------------
# rule to make .o from .f files
.f.o:
	$(FC) $(FFLAGS) $< -c
.f90.o:
	$(FC) $(FFLAGS) $< -c
#--------------------------------------
OBJ=\
gd2gc_cmdln.o\
gd2gc.o\
#--------------------------------------
$(main) : $(OBJ) $(main).o
	$(FC) $(FFLAGS) -o $(main) $(main).o $(OBJ)
#        rm *.o
#        rm *.mod
#--------------------------------------
