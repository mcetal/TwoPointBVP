FILES = ibvcond.f chtools.f discret.f int_density.f lu.f recur.f solve.f tree.f mesh.f misc.f bvp.f
OBJECTS = ibvcond.o chtools.o discret.o int_density.o lu.o recur.o solve.o tree.o mesh.o misc.o bvp.o
DIR = /home/math1/mkropins/Lib/bvp
BVP = libbvp.a

FFLAGS = -O

.f.o:	; f77 $(FFLAGS) -c $*.f

$(BVP): $(OBJECTS)
	ar r $(BVP) $(OBJECTS)
	ranlib $(BVP)
