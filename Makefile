# Compiler
CC = mpiicpx -mkl -O2 -xCORE-AVX512

# Flags
CFLAGS = -std=c++17 -I $(EIGEN3_INCLUDE_DIR)

# Source directories
SRC_DIR = src
INC_DIR = inc
OBJ_DIR = obj
BDIR = bin

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Object files
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Compile C++ source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
		$(CC) $(CFLAGS) -c $< -o $@ -I $(INC_DIR)

$(OBJ_DIR)/main.o:
		$(CC) $(CFLAGS) -c main.cpp -o $@ -I $(INC_DIR)

# Compile and link
main: $(OBJS) $(OBJ_DIR)/main.o
		$(CC) $(CFLAGS) $(OBJS) $(OBJ_DIR)/main.o -o $@ 

.PHONY: clean

clean:
		rm -f $(OBJS) $(BDIR)/main