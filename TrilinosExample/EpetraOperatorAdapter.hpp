/*
 * EpetraOperatorAdapter.hpp
 *
 *  Created on: Apr 10, 2012
 *      Author: uvilla
 */

#ifndef EPETRAOPERATORADAPTER_HPP_
#define EPETRAOPERATORADAPTER_HPP_

#include "EpetraVectorAdapter.hpp"
#include <Epetra_Operator.h>

class EpetraOperatorAdapter
{
public:
	EpetraOperatorAdapter(Epetra_Operator & op_) : op(op_) { };

	void Apply(const EpetraVectorAdapter & X, EpetraVectorAdapter & Y) const
	{
		int ierr(0);
		ierr = op.Apply(X.EpetraVector(), Y.EpetraVector() );

		assert(0 == ierr);
	}

	virtual ~EpetraOperatorAdapter(){ };

private:
	Epetra_Operator & op;
};

class EpetraPreconditionerAdapter
{
public:
	EpetraPreconditionerAdapter(Epetra_Operator & prec_) : prec(prec_) { };

	void Apply(const EpetraVectorAdapter & X, EpetraVectorAdapter & Y) const
	{
		int ierr(0);
		ierr = prec.ApplyInverse(X.EpetraVector(), Y.EpetraVector() );

		assert( 0 == ierr );
	}

	virtual ~EpetraPreconditionerAdapter(){ };

private:
	Epetra_Operator & prec;
};


#endif /* EPETRAOPERATORADAPTER_HPP_ */
