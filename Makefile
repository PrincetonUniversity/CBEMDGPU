CXX = icpc
PATHTOBOOST = /home/gkhoury/boost_1_52_0/
CFLAGS = -O2 -I $(PATHTOBOOST) 
OMPFLAGS = -openmp 
Default: MD
OMP = cellList.o integrator.o main.o nvt.o potential.o system.o utils.o

%.o : %.c
	$(CXX) $(CFLAGS) -c $<

cellList.o : cellList.cpp
	$(CXX) $(OMPFLAGS) $(CFLAGS) -c cellList.cpp
	
integrator.o : integrator.cpp
	$(CXX) $(OMPFLAGS) $(CFLAGS) -c integrator.cpp

main.o : main.cpp
	$(CXX) $(OMPFLAGS) $(CFLAGS) -c main.cpp
	
nvt.o : nvt.cpp
	$(CXX) $(OMPFLAGS) $(CFLAGS) -c nvt.cpp
	
potential.o : potential.cpp
	$(CXX) $(OMPFLAGS) $(CFLAGS) -c potential.cpp
	
system.o : system.cpp
	$(CXX) $(OMPFLAGS) $(CFLAGS) -c system.cpp

utils.o : utils.cpp
	$(CXX) $(OMPFLAGS) $(CFLAGS) -c utils.cpp

MD: $(OMP)
	$(CXX) $(OMPFLAGS) -o md $(CFLAGS) $(OMP) 

clean:
	$(RM) md
	$(RM) *.o
