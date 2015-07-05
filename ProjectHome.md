A templated parallel C++ implementation of the minres algorithm.

MINRES is a Krylov subspace iterative method for the solution of large sparse **symmetric** (possibly **indefinite**) **linear systems**.

The MINRES algorithm was first proposed by **C.C. Paige** and **M.A. Saunders** in 1975.

It is based on _Lanczos tridiagonalization_ and it allows for a _symmetric positive definite preconditioner_.

**Umberto Villa** (now at Lawrence Livermore National Laboratory) is the main author of the **tminres** project. **Santiago Akle** contributed to the project by providing the CUDA wrappers.

For more information regarding the MINRES algorithm, please visit Prof **M. Saunders** webpage at

http://www.stanford.edu/group/SOL/software/minres.html


---


The aim of **tminres** is to provide an _efficient_ and _portable_ implementation of the minres algorithm.

**tminres** supports user-defined **parallel** and **serial** data structure for the description of vectors,linear operators, and preconditioners.

The abstract classes ` Vector_trait ` and ` Operator_trait `  define the signatures of the methods that the user-defined classes should provide.

Users can decide whenever to use static (templates) or dynamic (derivated classes) polymorphism in order to provide such interfaces.

Two practical implementations of `Vector_trait` and `Operator_trait` are included in the code:
  * Serial code: A `VectorSimple` class implements the methods required by minres.
  * Parallel code: A `VectorEpetraAdapter` shows how to wrap parallel data structures in Trilinos (http://trilinos.sandia.gov/) to be used with tminres.