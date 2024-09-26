include ../../example.mk
include ../../common.mk

# Directory containing source files
SRC_DIR = src
BUILD_DIR = build

# List all source files (.cpp and .cu)
SRCS_CPP = $(wildcard $(SRC_DIR)/*.cpp)
SRCS_CU = $(wildcard $(SRC_DIR)/*.cu)

# Derive object files from source files, placing them in the build directory
OBJS_CPP = $(SRCS_CPP:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
OBJS_CU = $(SRCS_CU:$(SRC_DIR)/%.cu=$(BUILD_DIR)/%.o)

# Combine all object files
OBJS = $(OBJS_CPP) $(OBJS_CU)

# Add the TinyXML2 library to the LIBS variable
LIBS += -ltinyxml2

# Add C++17 standard flag
# CXX_STD = -std=c++17 

# Add flags
OPT += -DBOOST_ALLOW_DEPRECATED_HEADERS
# Define the target for the executable
sph_dlb: OPT := $(filter-out -DTEST_RUN,$(OPT))
sph_dlb: $(OBJS)
	$(CC) $(OPT) $(CXX_STD) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

# Define a target for the test executable if needed
sph_dlb_test: $(OBJS)
	$(CC) $(OPT) $(CXX_STD) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)



# Rule to compile .cu files to .o files in the build directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu | $(BUILD_DIR)
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)
	
# # Rule to compile .cpp files to .o files in the build directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CC) $(OPT) $(CXX_STD) -c -o $@ $< $(INCLUDE_PATH)

# # Rule to compile .cpp files to .o files in the build directory
# $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
# 	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)


# Ensure the build directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Default target
all: sph_dlb sph_dlb_test

# Target to run the sph_dlb executable
run: sph_dlb
	mpirun -np 4 ./sph_dlb

# Phony targets
.PHONY: clean all run

# Clean up build artifacts
clean:
	rm -f $(BUILD_DIR)/*.o *~ core sph_dlb sph_dlb_test 

# Clean up output files
cleanout:
	rm -f *.csv *.vtp *.pvtp *.txt
