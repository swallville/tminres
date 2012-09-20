/*
 * MMOperator.cu
 * Created on: Jun 26, 2012
 *  Author: Santiago Akle
 *
 */

#include "CusparseOperator.hpp"
#include "MMOperator.hpp"
#include "mmio.h"
#include "cusparse_v2.h"
#include <cuda_runtime.h>
#include <vector>
#include <stdio.h>

MMOperator::MMOperator(cusparseHandle_t handle_, FILE* mm_f): CusparseOperator(handle_)
{
   int m;
   MM_typecode mm_t;
   cudaError_t cudaStat;
   
   //Define the matrix descriptor
   cusparseStatus_t err = cusparseCreateMatDescr(&matDesc); 
   if(err!=CUSPARSE_STATUS_SUCCESS)
   {
      std::cerr << "Unable to allocate matrix descriptor\n";
      throw err;

   } 

   //Set the type to hermitian and the index base to one 
   cusparseSetMatType (matDesc, CUSPARSE_MATRIX_TYPE_SYMMETRIC);  
   cusparseSetMatFillMode(matDesc, CUSPARSE_FILL_MODE_UPPER);
   cusparseSetMatIndexBase (matDesc, CUSPARSE_INDEX_BASE_ONE) ;

   //Read the matrix description form the file
   int mm_io_err = mm_read_banner(mm_f,&mm_t);
   if(mm_io_err!=0)
   {
       std::cerr << "Error while reading matrix market file \n";
       throw mm_io_err;
   }
   if(mm_is_complex(mm_t))
   {

       std::cerr << "Complex matrices are not supported \n";
       throw TMINRES_CUBLAS_MM_UNSUPORTED;
   }
   if(mm_is_pattern(mm_t))
   {

       std::cerr << "Pattern matrices are not supported \n";
       throw TMINRES_CUBLAS_MM_UNSUPORTED;
   }
 
   if(!mm_is_sparse(mm_t))
   {
        std::cerr << "Dense matrices are not supported by this operator \n";
        throw TMINRES_CUBLAS_MM_UNSUPORTED;
   }
   if(!mm_is_matrix(mm_t))
   {
        std::cerr << "File must be a matrix to be used by this operator \n";
        throw TMINRES_CUBLAS_MM_UNSUPORTED; 
   }
   if(!mm_is_symmetric(mm_t))
   {
        std::cerr << "Matrix is not symmetric \n";
        throw TMINRES_CUBLAS_MM_UNSUPORTED; 
   }
   
   //This call will skip the comments and then read the size of the matrix
   int ret_code = mm_read_mtx_crd_size(mm_f,&m,&n,&nnz);
   if (ret_code !=0)
   {
     std::cerr << "Unable to read file size \n"; 
     throw TMINRES_CUBLAS_MM_UNSUPORTED;  
   }

   //Validate that the matrix is square  
   if(m!=n)
   {
     std::cerr << "Matrix is not square \n";
     throw TMINRES_CUBLAS_MM_UNSUPORTED;   
   } 

   //Matrix market sparse files are stored in COO format so
   //we need to transform the row_ix to csr
   std::vector<double> values(nnz);
   std::vector<int>    col_ix(nnz);
   std::vector<int>    row_ix(nnz);

   //Matrix market sparse file store the lower triangular part
   // in column major order. We require row major order to form 
   // the csr format, therefore we will interpret the column coordinates
   // as row coordinates and the row as column coordinates and assume that the
   // upper triangular section is stored. 
 
   int i;
   for(i = 0; i < nnz; ++i)
      fscanf(mm_f, "%d %d %lg\n", &col_ix[i], &row_ix[i], &values[i]);
 
   //Pointer to the row coordinates in the device
   int* cooRowIndA; 

   // Space for the row coordinates
   cudaStat = cudaMalloc((void**)&cooRowIndA,sizeof(int)*(nnz));
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }
   //Space for the row pointers
   cudaStat = cudaMalloc((void**)&csrRowPtrA,sizeof(int)*(n+1));
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }
   //Copy the rows and convert to CRS
   cudaStat = cudaMemcpy(cooRowIndA,&row_ix[0],sizeof(int)*row_ix.size(),cudaMemcpyHostToDevice);
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to copy data to device memory\n";
      throw cudaStat;
   }
   //Call the conversion routine
   cusparseStatus_t staus = cusparseXcoo2csr(handle,cooRowIndA,nnz,n,csrRowPtrA,CUSPARSE_INDEX_BASE_ONE);
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to covert to CSR format\n";
      throw cudaStat;
   }
    
   //Wait for the async call to finish
   cudaStat = cudaDeviceSynchronize();
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }
 
    // Space for the column coordinates
   cudaStat = cudaMalloc((void**)&csrColIndA,sizeof(int)*nnz);
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   } 
   //Load the column indices to the gpu
   cudaStat = cudaMemcpy(csrColIndA,&col_ix[0],sizeof(int)*col_ix.size(),cudaMemcpyHostToDevice);
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }
 
   //Free the space used for the row coordinates
   cudaStat = cudaFree(cooRowIndA);
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to free memory from the coo format conversion will proceed\n";
   }


   //Copy the values
   // Allocate the space for the vectors that define the matrix 
   cudaStat = cudaMalloc((void**)&csrValA,sizeof(double)*nnz);
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }

   //Load the matrix values to the gpu
   cudaStat = cudaMemcpy(csrValA,&values[0],sizeof(double)*values.size(),cudaMemcpyHostToDevice); 
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }
 
   
 
}

int MMOperator::get_n()
{
    return n;
}

int MMOperator::get_nnz()
{
    return nnz;
}
