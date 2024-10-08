# Compiler
CC = mpiicpc -mkl -O2 -xCORE-AVX512

# Flags
CFLAGS = -std=c++17 -I $(EIGEN3_INCLUDE_DIR) -w
LFLAGS = ./inc/pfapack/c_interface/libcpfapack.a  ./inc/pfapack/fortran/libpfapack.a

# Source directories
SRC_DIR = src
INC_DIR = inc
OBJ_DIR = obj
BDIR = bin

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Object files
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Compile and link
$(BDIR)/main: $(OBJS) $(OBJ_DIR)/main.o
		$(CC) $(CFLAGS) $(OBJS) $(OBJ_DIR)/main.o -o $@ $(LFLAGS)

# Compile C++ source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
		$(CC) $(CFLAGS) -c $< -o $@ -I $(INC_DIR) $(LFLAGS)

$(OBJ_DIR)/main.o:
		$(CC) $(CFLAGS) -c main.cpp -o $@ -I $(INC_DIR) $(LFLAGS)


.PHONY: clean

clean:
		rm -f $(OBJS) $(OBJ_DIR)/main.o $(BDIR)/main

