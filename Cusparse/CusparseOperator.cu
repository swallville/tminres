/*
 * CusparseOperator.cpp
 * Created on: Jun 26, 2012
 *  Author: Santiago Akle
 *
 */

#include "CusparseOperator.hpp"
#include "cusparse_v2.h"
#include <cuda_runtime.h>
#include <vector>
#include <time.h>
CusparseOperator::CusparseOperator(cusparseHandle_t handle_, std::vector<int> row_ptr, std::vector<int> col_ix, std::vector<double> vals): handle(handle_)
{
   cudaError_t cudaStat;
   //Set the number of non zeros
   nnz = vals.size(); 
   //Set the number of rows
   n   = row_ptr.size()-1;

   cusparseStatus_t err = cusparseCreateMatDescr(&matDesc); 
   if(err!=CUSPARSE_STATUS_SUCCESS)
   {
      std::cerr << "Unable to allocate matrix descriptor\n";
      throw err;

   }
    
   //Set the type to hermitian and the index base to one 
   cusparseSetMatType (matDesc, CUSPARSE_MATRIX_TYPE_GENERAL);  
   cusparseSetMatIndexBase (matDesc, CUSPARSE_INDEX_BASE_ONE) ;

   // Allocate the space for the vectors that define the matrix 
   cudaStat = cudaMalloc((void**)&csrValA,sizeof(double)*nnz);
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }
  // Space for the values
   cudaStat = cudaMalloc((void**)&csrColIndA,sizeof(int)*nnz);
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }
   //Space for the pointers to the row starts
   cudaStat = cudaMalloc((void**)&csrRowPtrA,sizeof(int)*(n+1));
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
 
   //Load the matrix values to the gpu
   cudaStat = cudaMemcpy(csrValA,&vals[0],sizeof(double)*vals.size(),cudaMemcpyHostToDevice); 
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }

   //Load the row start poitners to the gpu
   cudaStat = cudaMemcpy(csrRowPtrA,&row_ptr[0],sizeof(int)*row_ptr.size(),cudaMemcpyHostToDevice); 
   if(cudaStat != cudaSuccess)
   {
      std::cerr << "Unable to allocate device memory for operator\n";
      throw cudaStat;
   }
   
}

CusparseOperator::CusparseOperator(cusparseHandle_t handle_): handle(handle_)
{}

CusparseOperator::~CusparseOperator()
{
    cudaFree(csrValA);
    cudaFree(csrColIndA);
    cudaFree(csrRowPtrA);
    cusparseStatus_t cusparseStatus = cusparseDestroyMatDescr(matDesc);
    if(cusparseStatus != CUSPARSE_STATUS_SUCCESS)
    {
       std::cerr << "Unable to free descriptor\n";
       throw cusparseStatus;
    }
}


void CusparseOperator::Apply(const CusparseVector & x, CusparseVector & y) const
{
    const double alpha = 1;
    const double beta  = 0;
    cusparseStatus_t cusparseStatus = cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,n,n,nnz,&alpha,matDesc,csrValA,csrRowPtrA,csrColIndA,x.d_v,&beta,y.d_v);  
    if(cusparseStatus != CUSPARSE_STATUS_SUCCESS)
    {
       std::cerr << "Unable to execute matrix vector product, Error: "<< cusparseStatus<<"\n";
       throw cusparseStatus;
    }

}
