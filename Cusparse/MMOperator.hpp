/*
 * MMOpearator.hpp
 *
 *  Created on: Sep 15, 2012
 *      Author: Santiago Akle
 */

#ifndef MMOPERATOR_HPP_
#define MMOPERATOR_HPP_

#include "CusparseVector.hpp"
#include "CusparseOperator.hpp"
#include <cusparse_v2.h>

#define  TMINRES_CUBLAS_MM_UNSUPORTED 1;

class MMOperator: public CusparseOperator {
public:
	//! @class MMOperator
	/*!
	 * @brief  This class can be used to load matrix market symmetric matrices
     * and use cusparse to form products against it.
     * */

	//! Constructor
	/*!
	 * @param handle_ : Cusparse library handle 
     * @param mm_file : Pointer to the file in Matrix market format
     */	
    MMOperator(cusparseHandle_t handle_, FILE* mm_file);
    //! Returns the size of the loaded matrix
    int get_n();
    //! Returns the number of non zeros of the loaded matrix.
    int get_nnz();

};
#endif /* CUSPARSEOPERATOR_HPP_ */
