include ../../example.mk
include ../../common.mk

# Directory containing source files
SRC_DIR = src

# List all source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Derive object files from source files
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(SRC_DIR)/%.o)

# Add the TinyXML2 library to the LIBS variable
LIBS += -ltinyxml2

# Add C++17 standard flag
# CXX_STD = -std=c++17 -Wall 

# Add flags
OPT += -Wall

# Define the target for the executable
sph_dlb: OPT := $(filter-out -DTEST_RUN,$(OPT))
sph_dlb: $(OBJS)
	$(CC) $(OPT) $(CXX_STD) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

# Define a target for the test executable if needed
sph_dlb_test: $(OBJS)
	$(CC) $(OPT) $(CXX_STD) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

# Rule to compile .cpp files to .o files
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(OPT) $(CXX_STD) -c -o $@ $< $(INCLUDE_PATH)

# Default target
all: sph_dlb sph_dlb_test

# Target to run the sph_dlb executable
run: sph_dlb
	mpirun -np 4 ./sph_dlb

# Phony targets
.PHONY: clean all run

# Clean up build artifacts
clean:
	rm -f $(SRC_DIR)/*.o *~ core sph_dlb sph_dlb_test 

# Clean up output files
cleanout:
	rm -f *.csv	*.vtp *.pvtp *.txt
