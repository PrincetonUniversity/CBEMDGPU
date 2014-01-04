CXX = icpc
PATHTOBOOST = /home/gkhoury/boost_1_52_0/
CFLAGS = -O2 -I $(PATHTOBOOST) -g -Wall -pthread
OMPFLAGS = -openmp 

Default: MD

MD_DEPEND = cellList.o integrator.o nvt.o potential.o system.o utils.o 
OMP = main.o $(MD_DEPEND)
OMP_TESTS= unittests.o $(MD_DEPEND) gtest.a
OMP_TIMING = scaling_studies.o $(MD_DEPEND)
OMP_LMP = compare_lammps.o $(MD_DEPEND)

GTEST_DIR = /home/cdsilva/gtest-1.7.0
CPPFLAGS += -isystem $(GTEST_DIR)/include
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

unittests.o : unittests.cpp $(GTEST_HEADERS)
	$(CXX) $(OMPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c $<

%.o : %.c
	$(CXX) $(CFLAGS) -c $<

%.o : %.cpp
	$(CXX) $(OMPFLAGS) $(CFLAGS) -c $<

MD: $(OMP)
	$(CXX) $(OMPFLAGS) -o md $(CFLAGS) $^ 

TESTS: $(OMP_TESTS)
	$(CXX) $(OMPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o tests

TIMING: $(OMP_TIMING)
	$(CXX) $(OMPFLAGS) -o timing $(CFLAGS) $^

LMP_COMPARE: $(OMP_LMP)
	$(CXX) $(OMPFLAGS) -o lmp_compare $(CFLAGS) $^

clean:
	$(RM) md
	$(RM) tests
	$(RM) timing
	$(RM) lmp_compare
	$(RM) *.o
