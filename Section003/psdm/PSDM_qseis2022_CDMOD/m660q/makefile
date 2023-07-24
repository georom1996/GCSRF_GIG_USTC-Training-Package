########################################################################
#                                                                      #
# makefile script                                                      #
#                                                                      #
# Written by: Ling Chen,                                               #
# Last modified: Sat. 12, March 2005                                   #
#                                                                      #
########################################################################

#CC   = gcc  
CC   = gcc -O3 
#CC   = cc 
CC   = cc 
LIBS = -lm -Lgetpar/lib -lget
#F77  = g77 
#CC   = cc -g   
#F77  = f77 -vger 
#F77  = f77 
F77   = gfortran -O3 -ffixed-line-length-0 -m64
BIN  = ./


all: M660q_model

M660q_model: \
	M660q_model.o
	$(F77) -o $(BIN)M660q_model M660q_model.o

M660q_model.o: \
	M660q_model.for
	$(F77) -c M660q_model.for

clean: 
	rm -f \
	*.o

########################################################################
