CXX = mlicxx
CXXFLAGS = -fopenmp
SRC = main.cpp Rem.cpp SquareLattice.cpp
TARGET = run
MPIRUN = mpirun
NP = 2

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

run: $(TARGET)
	$(MPIRUN) -np $(NP) ./$(TARGET)

clean:
	rm -f $(TARGET)
