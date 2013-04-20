FSOURCES = filtering.f90 revDouble.f90
CSOURCES = fftw.c
CHEADERS = fftw.h
FOBJS = $(FSOURCES:.f90=.o)
COBJS =  $(CSOURCES:.c=.o)
LIBS  = -lfftw3 
FLAGS = -O2 -funroll-loops -fPIC
WARNINGS = -Wall 
SHOBJ = filtering.so

all: $(FOBJS) $(COBJS)
	gfortran --shared $(FOBJS) $(COBJS) -o $(SHOBJ) $(LIBS)
	cp  $(SHOBJ) ~/RLibs/

$(FOBJS): $(FSOURCES)
	gfortran $(WARNINGS) $(FLAGS) -c $(FSOURCES)

$(COBJS): $(CSOURCES) $(CHEADERS)
	gcc $(WARNINGS) $(FLAGS) -c $(CSOURCES)

clean:
	rm $(SHOBJ) *.o
