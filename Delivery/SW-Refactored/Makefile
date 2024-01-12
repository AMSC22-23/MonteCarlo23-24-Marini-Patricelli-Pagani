CXX = g++
CXXFLAGS = -std=c++17 -fopenmp
SRCS = main.cpp SquareLattice.cpp SwendsenWangParallel.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = SWprogram

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

