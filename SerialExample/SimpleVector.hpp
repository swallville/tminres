/*
 * SimpleVector.hpp
 *
 *  Created on: Apr 6, 2012
 *      Author: uvilla
 */

#ifndef SIMPLEVECTOR_HPP_
#define SIMPLEVECTOR_HPP_

#include <iostream>

class SimpleVector
{
public:
	//! @class SimpleVector
	/*!
	 * @brief Implementation of a dense serial vector according to Vector_traits
	 */

	//! Constructor
	/*!
	 * @param size int : the size of the vector
	 */
	SimpleVector(int size);
	//! Destructor
	virtual ~SimpleVector();

	//! Set all the entry of the Vector equal to val
	SimpleVector &  operator=(const double & val);
	//! Set the entry of the Vector equal to the entries in RHS
	SimpleVector & operator=(const SimpleVector & RHS);
	//! multiply THIS by a scalar value
	void Scale(const double & val);
	//! Create a new vector with the same structure of THIS. Values are not initialized.
	SimpleVector * Clone();

	//! Access entry i (non const version)
	double & operator[](const int i);
	//! Access entry i (const version)
	const double & operator[](const int i) const;
	//! Access entry i. if i < 0 return 0
	const double at(const int i) const;

	//! Fill the entries of the vector with random numbers. The vector is normalized with norm 1
	void Randomize(int seed);

	//! Print all the entries of the vector
	void Print(std::ostream & os);

	//! result = v1 + c2*v2
	friend void add(const SimpleVector & v1, const double & c2, const SimpleVector & v2, SimpleVector & result);
	//! result = c1*v1 + c2*v2
	friend void add(const double & c1, const SimpleVector & v1, const double & c2, const SimpleVector & v2, SimpleVector & result);
	//! result = alpha(v1 + v2)
	friend void add(const double & alpha, const SimpleVector & v1, const SimpleVector & v2, SimpleVector & result);
	//! result = v1 + v2 + v3
	friend void add(const SimpleVector & v1, const SimpleVector & v2, const SimpleVector & v3, SimpleVector & result);
	//! result = v1 - v2
	friend void subtract(const SimpleVector & v1, const SimpleVector & v2, SimpleVector & result);
	//! return the inner product of v1 and v2
	friend double InnerProduct(const SimpleVector & v1, const SimpleVector & v2);


private:
	double * vals;
	int size;
};



#endif /* SIMPLEVECTOR_HPP_ */
