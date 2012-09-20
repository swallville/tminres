CC = nvcc
C_FLAGS = -g -arch=sm_20 
L_FLAGS = -lcuda -lcusparse
I_DIR = -I ../ -I . 

all : cusparse_example3.exe mm_minres_benchmark.exe

cusparse_example3.exe: CusparseVector.o mmio.o CusparseOperator.o cusparse_example3.o 
	$(CC) $^ -o $@ -arch=sm_20 -lcusparse -lcuda

mm_minres_benchmark.exe: CusparseVector.o CusparseOperator.o mmio.o mm_minres_benchmark.o
	$(CC) $^ -o $@ -arch=sm_20 -lcusparse -lcuda

%.o : %.cu
	$(CC) -c $(C_FLAGS) $(I_DIR) $^ 

clean :
	rm *.o *.exe
