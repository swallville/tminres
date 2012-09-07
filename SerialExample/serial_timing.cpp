/*
 * serial_timing.cpp
 * Created: Jun 27, 2012
 * Author: Santiago Akle
 *
 */
 
/* This example generates a sequence of random problems 
 * which span a set of sizes. It solves each using Umberto Villas sequential
 * Each problem is timed and
 * the results are stored in results_out.txt
 */

#include <cstdlib>
#include "timing_minres.hpp"
#include "SimpleVector.hpp"
#include "SimpleOperator.hpp"
#include <fstream>
#include <cmath>
#include <ctime>
//Functor to generate the random data.
class rand_functor
{
public:
    rand_functor(int seed)
    {
        std::srand(seed);
    }
    rand_functor()
    {
    }
    double operator()()
    {
        return std::rand()/(double)RAND_MAX;
    }
};

//Fills a vector with random uniform numbers
void generate_random_values(std::vector<double> *d_a)
{

    std::generate(d_a->begin(),d_a->end(),rand_functor());
}

int main()
{
    int reps = 20;

    std::srand(0);
    //Initialize the cublas context

    //Instantiate the appropriate Cuda stuctures
    SimpleVector*   c_x;
    SimpleVector*   c_rhs;
    SimpleOperator* c_A;
    SimpleVector*   c_sol;
    
    //Pointers to hold the stl versions
    //SimpleVector*   s_x;
    //SimpleVector*   s_rhs;
    //SimpleOperator* s_A;
    
    std::vector<double> times;
    std::vector<double> errors;

    //File to log results
    std::ofstream log("serial_results.csv");
	
    int sizes[] = {5,10,20,40,80,160,320,640};
    for(int p = 0; p<8 ;p++)
    {
        for(int r = 0; r<reps; r++)
            {
              int n = sizes[p];
              
              c_A = new SimpleOperator(n);
              c_rhs = new SimpleVector(n);
              c_sol   = new SimpleVector(n);
              
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
              SimpleVector c_x(n);
              c_x = 0;
              
	          //(8) Solve the problem with minres
              int iterations;
              double resNorm;
              time_t tic = time(0);
	          MINRES<SimpleOperator,SimpleVector,SimpleOperator>(*c_A, c_x, *c_rhs, NULL, shift, max_iter, tol, show, &iterations,&resNorm);    
              time_t toc = time(0);

	          subtract(c_x, *c_sol, c_x);
	          double err2 = InnerProduct(c_x,c_x);

              delete c_sol;
              delete c_rhs;
              delete c_A;

              std::cout <<n<<","<<err2<<","<<iterations<<toc-tic<<","<<resNorm;
              log <<n<<","<<err2<<","<<iterations<<","<<resNorm<<","<<toc-tic;
              if(p!=7||r!=reps-1){log<<"\n";std::cout<<"\n";}
            }
	
    }

}

