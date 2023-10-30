# Makefile for DYNAMO command line subgraph search
# Patrick Mackey, 2018
SRC_DIR = src
OBJ_DIR = obj
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC))
CFLAGS = --std=c++14 -O2 -pthread -fopenmp
INCLUDES =
LDFLAGS = 
TARGET = graph_search

# Compiler (Must be g++ 4.9 or greater)
CXX = g++-9

.SUFFIXES:
.SUFFIXES: .o .cpp

all: $(TARGET)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CFLAGS)  $(INCLUDES) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $@

$(TARGET): $(OBJ)
	$(CXX) $(CFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

lint:
	clang-format $(SRC)

clean:
	rm -f $(OBJ) $(TARGET)

