########################################################################
#                                                                      #
# makefile script                                                      #
#                                                                      #
# Written by: Ling Chen,                                               #
# Last modified: Sep. 25, 2008                                         #
#                                                                      #
########################################################################

F77  = gfortran -ffixed-line-length-0 -m64 
BIN  = ./


all: pierc_new_n binr_vary_scan_n

pierc_new_n: \
	pierc_new_n.o utilities_nr.o
	$(F77) -o $(BIN)pierc_new_n pierc_new_n.o utilities_nr.o

pierc_new_n.o: \
	pierc_new_n.f
	$(F77) -c pierc_new_n.f

binr_vary_scan_n: \
	binr_vary_scan_n.o utilities_nr.o moveout_n.o utm_geo.o
	$(F77) -o $(BIN)binr_vary_scan_n binr_vary_scan_n.o utilities_nr.o moveout_n.o utm_geo.o

binr_vary_scan_n.o: \
	binr_vary_scan_n.f
	$(F77) -c binr_vary_scan_n.f

utilities_nr.o: \
	utilities_nr.f
	$(F77) -c utilities_nr.f

moveout_n.o: \
	moveout_n.f
	$(F77) -c moveout_n.f

utm_geo.o: \
	utm_geo.f
	$(F77) -c utm_geo.f

clean: 
	rm -f \
	*.o

########################################################################
