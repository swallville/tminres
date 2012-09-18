/*
 * CusparseOperator.cu
 *
 *  Created on: Aug 14, 2012
 *      Author: Santiago Akle
 */

#ifndef CUSPARSEOPERATOR_HPP_
#define CUSPARSEOPERATOR_HPP_

#include "CusparseVector.hpp"
#include <cusparse_v2.h>

class CusparseOperator {
public:
	//! @class CusparseOperator
	/*!
	 * @brief  This class represents a general sparse operator in CSR format
     * and implements the call to NVIDIA CUSPARSE to form products against it.
     * */

	//! Constructor
	/*!
	 * @param handle_ : Cusparse library handle 
     * @param row_ptr : Vector of length n+1 with the indices of the row starts plus nnz+1
     * @param col_ix  : Vector of length nnz with the column coordinates of the non zero values
     * @param vals    : Vector of length nnz with the values of the matrix
	 */
	
    CusparseOperator(cusparseHandle_t handle_, std::vector<int> row_ptr, std::vector<int> col_ix, std::vector<double> vals);

    ~CusparseOperator();

	//! y = A*x
	void Apply(const CusparseVector & X, CusparseVector & y) const;
    
protected:
    //!Constructor for inheriting classes
    /*
     * @param handle_ : Cusparse library handle
     */
    CusparseOperator(cusparseHandle_t handle_);

    //Cusparse session handle
    cusparseHandle_t handle;

    //Matrix descriptor
    cusparseMatDescr_t matDesc;

    //Size of the matrix
	int n;

    //Number of non zeros in the matrix
    int nnz;

    //Pointer to the vector containing the matrix values in the GPU
    double*  csrValA; 

    //Pointer to the vector containing the indices of the row starts in the GPU
    int*     csrRowPtrA;

    //Pointer to the vector containing the indices of the columns of the values in the GPU
    int*     csrColIndA;
};
#endif /* CUSPARSEOPERATOR_HPP_ */
