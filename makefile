#2013.10.07 G.S. Jung @LAMM
FFTW =/home/gsjung/applic/libs/fftw/lib 
VORO =/home/gsjung/applic/libs/voro/lib 
VORO_INC =/home/gsjung/applic/libs/voro/include
LMPS =/home/gsjung/applic/libs/lmp 
LMPS_INC=/home/gsjung/applic/include/lmp
REAX =/home/gsjung/Sources/lammps/lib/reax 


SRC = ./src
INC = ./src

FLAG = -g -d -I$(LMPS_INC)

#LDFLAG = -lm -L$(LMPS) -llmp -L$(FFTW) -lfftw3 -L$(VORO) -lvoro++ -L$(REAX) -lreax -L$(IFOR) -lifcore -lsvml -limf
LDFLAG = -lm -L$(LMPS) -llmp -L$(FFTW) -lfftw3 -L$(VORO) -lvoro++ -L$(REAX) -L$(IFOR) -lifcore -lsvml -limf -lcolvars
#LDFLAG = -lm -L$(LMPS) -L$(FFTW) -L$(VORO) -lvoro++ -L$(REAX) -L$(IFOR) -lifcore -lsvml -limf

TARGET = modeling.x 

CC = icc
CPP = icpc
MPICXX =mpicxx -I${VORO_INC}

OBJ_CPP = atominfo.o readfile.o lmpinput.o gensolid.o graphene.o fcc.o mympi.o grained.o

main.o : $(SRC)/main.cpp $(INC)/myheader.h $(OBJ_CPP)
	$(MPICXX) -c $(FLAG) $(SRC)/main.cpp

atominfo.o : $(SRC)/atominfo.cpp $(INC)/atominfo.h
	$(MPICXX) -c $(FLAG) $(SRC)/atominfo.cpp

readfile.o : $(SRC)/readfile.cpp $(INC)/readfile.h
	$(MPICXX) -c $(FLAG) $(SRC)/readfile.cpp

lmpinput.o : $(SRC)/lmpinput.cpp $(INC)/lmpinput.h
	$(MPICXX) -c $(FLAG) $(SRC)/lmpinput.cpp

gensolid.o : $(SRC)/gensolid.cpp $(INC)/gensolid.h
	$(MPICXX) -c $(FLAG) $(SRC)/gensolid.cpp

graphene.o : $(INC)/graphene.h $(SRC)/graphene.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/graphene.cpp

fcc.o : $(INC)/fcc.h $(SRC)/fcc.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/fcc.cpp

cellinfo.o : $(INC)/cellinfo.h $(SRC)/cellinfo.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/cellinfo.cpp

assemble.o : $(INC)/assemble.h $(SRC)/assemble.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/assemble.cpp

mympi.o : $(INC)/mympi.h $(SRC)/mympi.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/mympi.cpp

cell.o : $(INC)/cell.h $(SRC)/cell.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/cell.cpp

lmpcommand.o : $(INC)/lmpcommand.h $(SRC)/lmpcommand.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/lmpcommand.cpp

eam.o : $(INC)/eam.h $(SRC)/eam.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/eam.cpp

grained.o : $(INC)/grained.h $(SRC)/grained.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/grained.cpp

gslmp.o : $(INC)/gslmp.h $(SRC)/gslmp.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/gslmp.cpp

analysis.o : $(INC)/analysis.h $(SRC)/analysis.cpp
	$(MPICXX) -c $(FLAG) $(SRC)/analysis.cpp

bin  : $ main.o 
	$(MPICXX) -o $(TARGET) main.o $(OBJ_CPP) $(FLAGS) $(LDFLAG)

clean :
	rm -f  ${TARGET} *.o  *.dat log.* $(SRC)/*~ cfg/* tmp/*.tmp
