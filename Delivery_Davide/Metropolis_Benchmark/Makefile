# Compiler and flags
CXX = clang++
CXXFLAGS = -std=c++17 -fopenmp -I./

# Directories
SRC_DIR = .
OBJ_DIR = obj
BIN_DIR = bin

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp $(SRC_DIR)/**/*.cpp)

# Object files
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

# Executable
EXEC = $(BIN_DIR)/main

# Target for the executable
all: $(EXEC)

# Rule to compile object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to link the executable
$(EXEC): $(OBJS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Target to clean generated files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Target to run the executable
run: $(EXEC)
	$(EXEC)

.PHONY: all clean run
