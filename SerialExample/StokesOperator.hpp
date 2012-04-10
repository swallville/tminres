// tminres is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
//
// Authors:
// - Umberto Villa, Emory University - uvilla@emory.edu
// - Michael Saunders, Stanford University
// - Santiago Akle, Stanford University

/*!
@file
@author U. Villa - uvilla@emory.edu
@date 04/2012
*/

#ifndef STOKESOPERATOR_HPP_
#define STOKESOPERATOR_HPP_

#include "SimpleVector.hpp"

class StokesOperator {
public:
	//! @class StokesOperator
	/*!
	 * @brief A matrix free operator for the generalized Stokes problem
	 *
	 * This class implements a Markers and Cells (MAC) finite volume discretization of the generalized Stokes problem
	 * in the unit cube:
	 * $$ u - \Delta u + dp/dx = f_x $$
	 * $$ v - \Delta v + dp/dy = f_y $$
	 * $$ du/dx + dv/dy        = 0   $$
	 *
	 * where [u,v] is the velocity field and p the pressure fields.
	 */

	//! Constructor
	/*!
	 * @param n_ int : number of cell in the x direction. The number of cell in the y direction is also n_.
	 */
	StokesOperator(int n_ );

	//! Do nothing distructor.
	virtual ~StokesOperator();

	//! Y = A*X where A is the Finite Volumes matrix of the discrete system.
	void Apply(const SimpleVector & X, SimpleVector & Y) const;

private:

	inline int getUIndex(int i, int j) const
	{
		if(i == -1)
			i = n;

		if(i==n+1)
			i = 0;

		if(j== -1)
			j = n-1;

		if(j==n)
			j = 0;

		return offsetU + i + (n+1)*j;
	}

	inline int getVIndex(int i, int j) const
	{
		if(i == -1)
			i = n-1;

		if(i==n)
			i = 0;

		if(j== -1)
			j = n;

		if(j==n+1)
			j = 0;

		return offsetV + i + n*j;
	}

	inline int getPIndex(int i, int j) const
	{
		if(i == -1 || i == n || j == -1 || j == n)
			return -1;

		return offsetP + i + n*j;
	}


	int n;
	double h;
	int offsetU;
	int offsetV;
	int offsetP;
};

class StokesDiagonalScaling
{
public:
	//! @class StokesDiagonalScaling
	/*!
	 * @brief A simple Diagonal preconditioner for the generalized Stokes Problem
	 */

	//! Constructor
	/*!
	 * @param n_ int : number of cell in the x direction. The number of cell in the y direction is also n_.
	 */
	StokesDiagonalScaling(int n_);

	virtual ~StokesDiagonalScaling();

	//! Y = M\X where M is the diagonal of the matrix in StokesOperator.
	void Apply(const SimpleVector & X, SimpleVector & Y) const;

private:
	int n;
	double h;
	int offsetU;
	int offsetV;
	int offsetP;

};

#endif /* STOKESOPERATOR_HPP_ */
