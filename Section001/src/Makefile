#created by Joe zhou in guangzhou 20200319, contact with me on zhangzhou1996@gmail.com
CCC=g++
FFF=gfortran
F90=gfortran
SUBS=sac_lpz.o
CFLAG=-std=c99
# CFLAG=
SUMV=subs.o saclib.o
all:  modules1 modules2 SX_Deconv sacStack GC_SRF GC_SRF_par GC_SRF_post spmove.joe sacsnr
#Compile modules.
modules1 :
	$(F90) -O4 -c mod_sac_io.f90
modules2 : subs.f saclib.f
	$(FFF) -c subs.f saclib.f
SX_Deconv: 
	$(CCC) $(CFLAG) -c sac_lpz.c 
	$(CCC)  -g -c $@.cpp 
	$(CCC) $(CFLAG) -o $@ $@.o sac_lpz.o
	mv $@ ../bin/
sacStack:
	gcc $(CFLAG) -c sac_lpz.c sac_lpz.c fft.c Complex.c
	gcc $(CFLAG) -g -c $@.c
	gcc $(CFLAG) -o $@ $@.o sac_lpz.o fft.o Complex.o -lm
	mv $@ ../bin/
GC_SRF: 
	$(CCC) $(CFLAG) -c sac_lpz.c
	$(CCC)  -g -c $@.cpp
	$(CCC) $(CFLAG) -o $@ $@.o sac_lpz.o
	mv $@ ../bin/
GC_SRF_par:
	mpic++ $(CFLAG) -c sacio.c
	mpic++ $(CFLAG) -g -o $@ $@.cpp $(SUBS)
	mv $@ ../bin/	
GC_SRF_post:
	mpic++ $(CFLAG) -c sacio.c
	mpic++ -g -c $@.cpp
	mpic++ $(CFLAG) -o $@ $@.o sacio.o		
	mv $@ ../bin/
spmove.joe: modules2
	$(FFF) -c spmove.joe.f
	$(FFF) -o $@ $@.o subs.o saclib.o
	mv $@ ../bin/
sacsnr : modules1
	$(F90) -O4 sacsnr.f90 -o sacsnr ./mod_sac_io.o
	mv sacsnr ../bin/

clean:
# rm GC_SRF GC_SRF_par GC_SRF_post SX_Deconv spmove.joe
	rm -f *.o ../bin/*
