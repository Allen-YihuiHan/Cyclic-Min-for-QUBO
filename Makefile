# NVIDIA C compiler and flags.
CC = gcc
NVCC = nvcc

LDFLAGS	= -L$(CUDA_ROOT)/lib64 -lcudart

# Build up flags and targets
NVCCFLAGS = -O3 -Xcompiler -march=native

NVCCFLAGS += -DWITH_CUBLAS -I $(CUDA_ROOT)/include
LDFLAGS += -lcublas -L $(CUDA_ROOT)/lib64 -Xcompiler "-Wl,-rpath,$(CUDA_ROOT)/lib64"

EXECUTABLE = host

CU_SOURCES = host.cu kernel.cu 
C_SOURCES = binary_search_test.c genetic_algorithm.c

CU_OBJECTS = $(CU_SOURCES:.cu=.o)
C_OBJECTS = $(C_SOURCES:.c=.o)

all: $(EXECUTABLE)

$(CU_OBJECTS): %.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(C_OBJECTS): %.o: %.c
	$(CC) -c $< -o $@

$(EXECUTABLE): $(CU_OBJECTS) $(C_OBJECTS)
	$(NVCC) $^ -o $@ $(LDFLAGS)

clean:
	rm -f *.o