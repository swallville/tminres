/*
 * main.cpp
 *
 *  Created on: Apr 7, 2012
 *      Author: uvilla
 */

/*
 * A simple example of minres usage with diagonal preconditioner.
 * We consider a Markers and Cells Finite Volume discretization of a 2d Stokes problem:
 * u - Laplacian u + grad p = f in [0,1]x[0,1]
 * div u = 0 in [0,1]x[0,1]
 * where u is the velocity field and p the pressure field.
 *
 */

#include <pminres.hpp>
#include "SimpleVector.hpp"
#include "StokesOperator.hpp"
#include <cmath>

int main()
{
	//(1) Discretization parameters and problem size
	int n(4);					//number of cell for each edge
	int dim(n*n + 2*(n+1)*n);	//total number of unknows (n^2 pressure and (n+1)n for each velocity component)
	//(2) Define the linear operator "op" we want to solve.
	StokesOperator op(n);
	//(3) Define the exact solution (at random)
	SimpleVector sol(dim);
	sol.Randomize( 1 );

	//(4) Define the "rhs" as "rhs = op*sol"
	SimpleVector rhs(dim);
	op.Apply(sol, rhs);
	double rhsNorm( sqrt(InnerProduct(rhs,rhs)) );
	std::cout << "|| rhs || = " << rhsNorm << "\n";

	//(5) Define the preconditioner
	StokesDiagonalScaling prec(n);

	//(6) Use an identically zero initial guess
	SimpleVector x(dim);
	x = 0;

	//(7) Set the minres parameters
	double shift(0);
	int max_iter(100);
	double tol(1e-6);
	bool show(true);

	//(8) Solve the problem with minres
	MINRES(op, x, rhs, &prec, shift, max_iter, tol, show);

	//(9) Compute the error || x_ex - x_minres ||_2
	subtract(x, sol, x);
	double err2 = InnerProduct(x,x);

	std::cout<< "|| x_ex - x_n || = " << sqrt(err2) << "\n";

	return 0;
}
