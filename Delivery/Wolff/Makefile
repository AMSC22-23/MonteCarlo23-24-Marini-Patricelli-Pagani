CXX = g++
CXXFLAGS = -std=c++17 -fopenmp
TARGET = Wprogram

SRCS = main.cpp SquareLattice.cpp Wolff.cpp
OBJS = $(SRCS:.cpp=.o)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)

