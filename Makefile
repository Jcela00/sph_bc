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
LIBS_SELECT += -ltinyxml2

# Add flags if we're using NVCC
ifdef CUDA_CC
    ifeq ($(CUDA_CC), nvcc -ccbin=mpic++)
        CUDA_OPTIONS += -Xcompiler -Wno-attributes -diag-suppress=549 -diag-suppress=68 -DBOOST_ALLOW_DEPRECATED_HEADERS
    endif
endif

CUDA_OPTIONS := $(filter-out -DTEST_RUN, $(CUDA_OPTIONS))

# Rule to compile .cu files to .o files in the build directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu | $(BUILD_DIR)
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

# Rule to compile .cpp files to .o files in the build directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CC) $(OPT) $(CXX_STD) -c -o $@ $< $(INCLUDE_PATH)

# Build the executable
sph: $(OBJS)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

# Ensure the build directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Default target
all: sph sph_test

# Target to run the sph executable
run: sph 
	./sph XML/example/Cavity.xml

# Phony targets
.PHONY: clean all run

# Clean up build artifacts
clean:
	rm -f $(BUILD_DIR)/*.o *~ core sph sph_test 

# Clean up output files
cleanout:
	rm -f *.csv *.vtp *.pvtp *.txt
print_cuda_cc:
	@echo "CUDA_CC is set to: $(CUDA_CC)"