# Compiler
CC = mpiicpc -mkl -O2 -xCORE-AVX512

# Flags
CFLAGS = -std=c++17 -I /DATA.c9/wanzhouquan/V-PQMC/cPQMC/eigen-3.4.0

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

# Compile and link
main: $(OBJS)
		$(CC) $(CFLAGS) $(OBJS) -o $@ 

.PHONY: clean

clean:
		rm -f $(OBJS) $(BDIR)/main