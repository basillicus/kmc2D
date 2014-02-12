#DEBUG   = -nbs -C -d2
FFLAGC  =   -O2
#FFLAGC  =  -O -C
#FFLAGC  =  -g
opt     = 
profile = 
obj     = kmc
ildir   =  
#LIBS    = -L/Applications/Absoft11.0/lib64/ -lU77
input   = 

CFLAGC  = -O

FCOMPL  = gfortran 
#FCOMPL  = ifort -debug -check -traceback -fpe0 -warn -Vaxlib $(DEBUG)
#FCOMPL  = ifort  -Vaxlib $(DEBUG)

# list of other directories for source files
.PREFIXES: .
.SUFFIXES:
.SUFFIXES: .f90 .c .s .o .fil 

.f90.o:
	$(FCOMPL) -c $(ildir) $(FFLAGC) $(opt) $(profile) $<

.f.fil:
	$(FCOMPL) -il $(FFLAGC) $<

.s.o:
	as $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) $(profile) $<

OBJECTS = change.o moves.o tools_strings.o param.o info.o interactions.o  lattice.o  occupied.o draw.o make_move.o tools.o main.o random.o hat.o objects.o 

INLINE  =

APPLIC: $(INLINE) $(OBJECTS) $(OBJECTS1)
	$(FCOMPL) $(FFLAGC) $(profile) -o $(obj) $(OBJECTS) $(LIBS)

test:
	@echo START TEST ON $(input) , opt = $(opt)
	@echo start test on $(input) , opt = $(opt) >> TIME.LOG
	@date >> TIME.LOG
	@( time $(obj) < $(input) > $(input).out ) 2>> TIME.LOG
	@echo - - - - - - - - - - - >> TIME.LOG

clean:
	@rm -f $(INLINE) $(OBJECTS) $(OBJECTS1) $(OBJ_mod) *.mod

# include file dependencies
param.o :: tools_strings.o
info.o :: param.o lattice.o
tools.o :: info.o  
occupied.o :: info.o moves.o
main.o :: occupied.o change.o 
make_move.o :: moves.o info.o change.o 

draw.o :: info.o 

hat.o :: build.hat



