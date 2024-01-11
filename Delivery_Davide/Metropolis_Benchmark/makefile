CXX = clang++
CXXFLAGS = -std=c++17 -fopenmp -MMD -MP
SRCDIR = /Users/davidepatricelli/Desktop/src
METROPOLIS_DIR = $(SRCDIR)/Metropolis
LATTICE_DIR = $(SRCDIR)/Lattice
INCLUDES = -I$(METROPOLIS_DIR) -I$(LATTICE_DIR)

SOURCES = $(wildcard $(METROPOLIS_DIR)/*.cpp) $(wildcard $(LATTICE_DIR)/*.cpp) $(SRCDIR)/main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
DEPENDS = $(OBJECTS:.o=.d)
EXECUTABLE = program

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

-include $(DEPENDS)

clean:
	rm -f $(OBJECTS) $(DEPENDS) $(EXECUTABLE)
