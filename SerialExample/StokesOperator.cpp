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

#include "StokesOperator.hpp"
#include <cassert>

StokesOperator::StokesOperator(int n_):
	n(n_),
	h(1./static_cast<double>(n_)),
	offsetU(0),
	offsetV( (n+1)*n ),
	offsetP( (n+1)*n + n*(n+1) )
{
	assert(n > 0);
}

StokesOperator::~StokesOperator()
{
	// TODO Auto-generated destructor stub
}

void StokesOperator::Apply(const SimpleVector & X, SimpleVector & Y) const
{
	Y = 0.;
	for(int i(0); i<n+1; ++i)
		for(int j(0); j<n+1; ++j)
		{
			double uij ( X[getUIndex(i  ,j  ) ] );
			double uimj( X[getUIndex(i-1,j  ) ] );
			double uipj( X[getUIndex(i+1,j  ) ] );
			double uijm( X[getUIndex(i  ,j-1) ] );
			double uijp( X[getUIndex(i  ,j+1) ] );

			double vij ( X[getVIndex(i  ,j  ) ] );
			double vimj( X[getVIndex(i-1,j  ) ] );
			double vipj( X[getVIndex(i+1,j  ) ] );
			double vijm( X[getVIndex(i  ,j-1) ] );
			double vijp( X[getVIndex(i  ,j+1) ] );

			double pij ( X.at(getPIndex(i  ,j  ) ) );
			double pimj( X.at(getPIndex(i-1,j  ) ) );
			double pijm( X.at(getPIndex(i  ,j-1) ) );

			double laplacianu( 4.*uij - uimj - uipj - uijm - uijp );
			laplacianu *= (1/h/h);
			double laplacianv( 4.*vij - vimj - vipj - vijm - vijp );
			laplacianv *= (1/h/h);

			double dxp( pij - pimj);
			dxp *= 1./h;
			double dyp( pijm - pij);
			dyp *= 1./h;

			double div( uipj - uij + vij - vijp);
			div *= -1./h;

			if ( j != n )
				Y[getUIndex(i,j)] += laplacianu + dxp + uij;

			if ( i !=n )
				Y[getVIndex(i,j)] += laplacianv + dyp + vij;

			if( i != n && j != n )
				Y[getPIndex(i,j)] += div;


		}
}


StokesDiagonalScaling::StokesDiagonalScaling(int n_):
	n(n_),
	h(1./static_cast<double>(n_)),
	offsetU(0),
	offsetV( (n+1)*n ),
	offsetP( (n+1)*n + n*(n+1) )
{
	assert(n > 0);
}


StokesDiagonalScaling::~StokesDiagonalScaling()
{

}

void StokesDiagonalScaling::Apply(const SimpleVector & X, SimpleVector & Y) const
{
	int i(0);

	double diagValue(4./h/h + 1.);
	for( ; i < offsetP; ++i)
		Y[i] = X[i]/diagValue;

	for( ; i < offsetP + n*n; ++i)
		Y[i] = X[i];
}
