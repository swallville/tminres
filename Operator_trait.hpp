/*
 * Operator_trait.hpp
 *
 *  Created on: Apr 10, 2012
 *      Author: uvilla
 */

#ifndef OPERATOR_TRAIT_HPP_
#define OPERATOR_TRAIT_HPP_

#include "Vector_trait.hpp"

class Operator_trait
{
	//! @class Operator_trait
	/*!
	 * @brief This class defines the interface of a linear operator to be used in MINRES
	 */
public:
	//! Y = A * X
	void Apply(const Vector_trait & X, Vector_trait & Y) const = 0;
};

class Preconditioner_trait
{
	//! @class Preconditioner_trait
	/*!
	 * @brief This class defines the interface of a linear operator to be used as a preconditioner in MINRES
	 */
public:
	//! Y = M \ X;
	void Apply(const Vector_trait & X, Vector_trait & Y) const = 0;
};


#endif /* OPERATOR_TRAIT_HPP_ */
