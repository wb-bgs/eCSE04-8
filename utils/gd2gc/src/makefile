# decide which compiler to use the library depending on the architecture
FC = ftn
FFLAGS = -O3

#--------------------------------------
# rule to make .o from .f files
%.o : %.f
	$(FC) $(FFLAGS) $< -c
%.o : %.f90
	$(FC) $(FFLAGS) $< -c
#--------------------------------------
OBJ=\
gd2gc.o\
gd2gc_cmdln.o\
#--------------------------------------
gd2gc_cmdln : $(OBJ) gd2gc_cmdln.o
	$(FC) $(FFLAGS) -o gd2gc_cmdln $(OBJ)
#        rm *.o
#        rm *.mod
#--------------------------------------
