EXEC   = laces_create_sed

OPTIMIZE  =  -O2

OBJS   = main.o constants.o routines.o bpass_models.o igm-absorption.o madau_absorption.o pass_bands.o nebular_line_emission.o

CC     = g++

INCL   = constants.h routines.hpp nebular_line_emission.h bpass_models.h pass_bands.h igm-absorption.h madau_absorption.h

LIBS   = -lm -lgsl -lgslcblas 

CFLAGS   = $(OPTIMIZE)
CXXFLAGS = $(OPTIMIZE)


$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LIBS) -o $(EXEC)   

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

