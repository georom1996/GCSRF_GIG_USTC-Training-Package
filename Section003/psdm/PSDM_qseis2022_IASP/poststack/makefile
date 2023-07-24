FC = gfortran -ffixed-line-length-0
FFLAGS =  -O3 

OBJ = hdpmig.o fourt.o pdpm_sub.o phs_scrn.o trigenh.o trisolv.o wide_cpns.o xtaper.o 

fd: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o hdpmig.x

def:
	$(FC) $(FFLAGS) -c *.f
.f.o: $*.f
	$(FC) $(FFLAGS) -c $*.f

clean:
	/bin/rm -f *.o core
