/*
 * cuda_timing.cpp
 * Created: Jun 27, 2012
 * Author: Santiago Akle
 *
 */
 
/* This example generates a sequence of random problems 
 * which span a set of sizes. It solves each using Umberto Villas sequential
 * operator and Akles Cuda extension. Each problem is timed and
 * the results are stored in results_out.txt
 */

#include "timing_minres.hpp"
#include "CudaVector.hpp"
#include "CudaOperator.hpp"
#include <fstream>
#include <cmath>
#include <ctime>
//Functor to generate the random data.
class rand_functor
{
public:
    rand_functor(int seed)
    {
        srand(seed);
    }
    rand_functor()
    {
    }
    double operator()()
    {
        return rand()/(double)RAND_MAX;
    }
};

//Fills a vector with random uniform numbers
void generate_random_values(std::vector<double> *d_a)
{

    std::generate(d_a->begin(),d_a->end(),rand_functor());
}

//Generates a random linear operator, a right hand side of unit norm and a solution
void generate_random_problem(cublasHandle_t handle, int n, CudaVector* c_x_o, CudaVector* c_rhs, CudaOperator* c_A)
//void generate_random_problem(cublasHandle_t handle, int n, CudaVector* c_x_o, CudaVector* c_rhs, CudaOperator* c_A, SimpleVector* s_x, SimpleVector* s_rhs, SimpleOperator* s_A)
{
   
    std::vector<double> operator_dat(n*n);
    std::vector<double> x_dat(n); 
    
    generate_random_values(&operator_dat);
    generate_random_values(&x_dat);
    
    //Generate the CUDA objects
    c_A    = new CudaOperator(handle,n,operator_dat);
    c_x_o  = new CudaVector(handle,n);
    c_rhs  = new CudaVector(handle,n);
    *c_x_o = x_dat;
    //Generate the rhs 
    c_A->Apply(*c_x_o,*c_rhs);
    double rhs_norm = InnerProduct(*c_rhs,*c_rhs); 
    c_rhs->Scale(1./rhs_norm);
    c_x_o->Scale(1./rhs_norm); 

    /*
    //Generate the serial objects
    s_A    = new SimpleOperator(operator_dat); 
    s_x    = new SimpleVector(n); 
    s_rhs  = new SimpleVector(n);
    *s_x   = x_dat;
    s_A->Apply(*s_x,*s_rhs);
    s_rhs->Scale(1./rhs_norm);
    s_x->Scale(1./rhs_norm);
    */
    

}

int main()
{
    int reps = 20;

    srand(0);
    //Initialize the cublas context
    cublasHandle_t handle;    
    cublasStatus_t status = cublasCreate(&handle);
    if(status != CUBLAS_STATUS_SUCCESS)
    {
        std::cerr << "Unable to initialize cublas! \n";
        throw status;
    }
    std::cout << "Initialized cublas\n";

    //Instantiate the appropriate Cuda stuctures
    CudaVector*   c_x;
    CudaVector*   c_rhs;
    CudaOperator* c_A;
    CudaVector*   c_sol;
    
    //Pointers to hold the stl versions
    //SimpleVector*   s_x;
    //SimpleVector*   s_rhs;
    //SimpleOperator* s_A;
    
    std::vector<double> times;
    std::vector<double> errors;

    //File to log results
    std::ofstream log("cuda_results.csv");
	
    int sizes[] = {5,10,20,40,80,160,320,640};
    for(int p = 0; p<8 ;p++)
    {
        for(int r = 0; r<reps; r++)
            {
              int n = sizes[p];
              
              c_A = new CudaOperator(handle,n);
              c_rhs = new CudaVector(handle,n);
              c_sol   = new CudaVector(handle,n);
              
              std::vector<double> rand_vals(n);
              generate_random_values(&rand_vals);
              *c_sol = rand_vals;
              
              std::vector<double> rand_A_vals(n*n);
              generate_random_values(&rand_A_vals);
              *(c_A)=rand_A_vals;
              c_A->Apply(*c_sol,*c_rhs);
              
              double norm = InnerProduct(*c_rhs,*c_rhs);
              norm = sqrt(norm);
              c_sol->Scale(1.0/norm);
              c_rhs->Scale(1.0/norm);
              
              //Call minres
              double shift(0);
	          int max_iter(2*n);
	          double tol(1e-6);
	          bool show(false);
   
              //Store the solution 
              CudaVector c_x(handle,n);
              c_x = 0;
              
	          //(8) Solve the problem with minres
              int iterations;
              double resNorm;
              time_t tic = time(0);
	          MINRES<CudaOperator,CudaVector,CudaOperator>(*c_A, c_x, *c_rhs, NULL, shift, max_iter, tol, show, &iterations,&resNorm);    
              time_t toc = time(0);

	          subtract(c_x, *c_sol, c_x);
	          double err2 = InnerProduct(c_x,c_x);

              delete c_sol;
              delete c_rhs;
              delete c_A;

              log <<n<<","<<err2<<","<<iterations<<","<<resNorm<<","<<toc-tic;
              if(p!=7||r!=reps-1){log<<"\n";std::cout<<"\n";}
            }
	
    }

}


