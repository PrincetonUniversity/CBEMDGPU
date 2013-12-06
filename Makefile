CXX = g++
PATHTOBOOST = /scratch/pvfs2/gkhoury/boost_1_52_0/
CFLAGS = -O2 -pedantic -I $(PATHTOBOOST) 
Default: MD
OMP = cellList.o integrator.o main.o nvt.o potential.o system.o utils.o

%.o : %.c
	$(CXX) $(CFLAGS) -c $<

cellList.o : cellList.cpp
	$(CXX) $(CFLAGS) -c cellList.cpp
	
integrator.o : integrator.cpp
	$(CXX) $(CFLAGS) -c integrator.cpp

main.o : main.cpp
	$(CXX) $(CFLAGS) -c main.cpp
	
nvt.o : nvt.cpp
	$(CXX) $(CFLAGS) -c nvt.cpp
	
potential.o : potential.cpp
	$(CXX) $(CFLAGS) -c potential.cpp
	
system.o : system.cpp
	$(CXX) $(CFLAGS) -c system.cpp

utils.o : utils.cpp
	$(CXX) $(CFLAGS) -c utils.cpp

MD: $(OMP)
	$(CXX) -fopenmp -lgomp -o md $(CFLAGS) $(OMP) 

clean:
	$(RM) md
	$(RM) *.o
