/* This program generates p random problems for each
   matrix market file in the ./MM directory.
   It solves and logs the name of the file, the problem instance
   the iteration count the residual norm and the error.
*/

/*
 * cusparse_example3.cu
 * Created: September 16, 2012
 * Author: Santiago Akle
 *
 */
 
/* This example loads the matrix in matrix market format and solves 
 * using a random rhs.
 */
#include "timing_minres.hpp"
#include "CusparseVector.hpp"
#include "CusparseOperator.hpp"
#include <cmath>
#include <stdio.h>
#include <dirent.h>
#include <string>
#include <fstream>
#include <ctime>
#include <numeric>

#define REPS 20

class rand_functor
{
public:
    //Generates uniform random variables in the -1,1 range
    rand_functor(int seed)
    {
        srand(seed);
        std_f = 1;
    }

    //Generates uniformly distributed variables in a symmetric range such that 
    // the resulting standard deviation is std_
    rand_functor(int seed,double std_)
    {
        srand(seed);
        std_f = sqrt(3./7.)*std_;
    }
    double operator()()
    {    
        return std_f*((2*rand()/(double)RAND_MAX) -1);
    }
private:
    //Standard deviation factor
    double std_f;
};

void generate_random_values(std::vector<double> &d_a)
{

    std::generate(d_a.begin(),d_a.end(),rand_functor(0));
}

void generate_random_normalized_values(std::vector<double> &d_a)
{
    std::generate(d_a.begin(),d_a.end(),rand_functor(0,1./sqrt(d_a.size())));
}

void run_test_for_matrix(cusparseHandle_t handle, std::string filename, std::ofstream &results_file)
{

    std::cout << "Will load "<<filename<<"\n";
    FILE* mm_f = fopen(filename.c_str(),"r");
    CusparseOperator mm_O = CusparseOperator::LoadMMFile(handle,mm_f);
    int n = mm_O.get_n();
    int nnz = mm_O.get_nnz();

    fclose(mm_f);

    std::vector<double> dat(n);
    std::vector<double> errors;
    std::vector<double> residuals;
    std::vector<int>    iterations;
    std::vector<int>    times;

    int i;
    for(i=0 ; i<REPS ; i++)
    {
        //Generate a random rhs
        generate_random_normalized_values(dat);
        CusparseVector x_o(handle,n);
        CusparseVector rhs(handle,n);
        x_o = dat;
        mm_O.Apply(x_o,rhs);
        //The problem has been generated
        
        //Call timing minres
        int    iters; 
        double resnorm;
        double err2;

        //Set minres parameters
        double shift(0);
        int max_iter(10000);
        double tol(1e-6);
        bool show(false);
        
        CusparseVector x(handle,n);
        x = 0;
        
        clock_t tic = clock();
        //(8) Solve the problem with minres
        MINRES<CusparseOperator,CusparseVector,CusparseOperator>(mm_O, x, rhs, NULL, shift, max_iter, tol, show,&iters,&resnorm);
        clock_t toc = clock();

        //Calculate the error norm
        subtract(x, x_o, x);
	    err2 = InnerProduct(x,x);
        

        //XXX:Do we have to delete the vectors or will the change of context call the destructor?

        //Save the values
        std::cout << filename << "," << i << "," << iters << "," << resnorm << "," << sqrt(err2) << "," << (int)(toc-tic) << std::endl;
        errors.push_back(sqrt(err2));
        iterations.push_back(iters);
        residuals.push_back(resnorm);
        times.push_back(toc-tic);

    }
        double avg_error = std::accumulate(errors.begin(),errors.end(),0.0);
        avg_error /= REPS;
        double avg_resnorm = std::accumulate(residuals.begin(),residuals.end(),0.0);
        avg_resnorm /= REPS;
        double avg_iterations = std::accumulate(iterations.begin(),iterations.end(),0.0);
        avg_iterations /= REPS;
        double avg_time = std::accumulate(times.begin(),times.end(),0.0);
        avg_time /= REPS;

        results_file<< filename << "," << n <<"," << nnz <<"," << avg_iterations << "," << avg_resnorm << "," << avg_error << "," << avg_time << std::endl;
        results_file.flush();

}

int main()
{
  char* dirname = "../MatrixMarket/";
  DIR* dir = opendir(dirname);
  DIR* dir_c;
  struct dirent *ent;
  if(dir == NULL)
  {
      std::cerr << "Unable to open dir\n";
      exit(-1);
  }
  std::vector<std::string> files;
  std::vector<std::string> dirs;
  while((ent = readdir(dir)) != NULL)
  {
      if(ent->d_type == DT_DIR)
      { 
          std::string sub_dirname = std::string(ent->d_name);
          if(sub_dirname!=".")
              if(sub_dirname!="..")
                dirs.push_back(std::string(dirname)+sub_dirname+"/");
      }
  }
  for(int j=0;j<dirs.size();j++)
  {
     dir_c = opendir(dirs[j].c_str()); 
     while((ent = readdir(dir_c))!=NULL)
     {
         if(ent->d_type !=DT_DIR)
         {
            //Check the termination 
             std::string file_name  =std::string(ent->d_name);
            int file_len = file_name.size();
            //Select the files that end in mtx
            if(file_name.compare(file_len-4,4,".mtx")==0){
                std::cout << "Will benchmark against " << file_name << std::endl;
                files.push_back(dirs[j]+file_name);
            }
         }
     }
     closedir(dir_c);
  }
  std::cout << "List of MM files:\n";
  for(int j = 0;j<files.size();j++)
      std::cout << files[j]<< std::endl;
  closedir(dir);

  //Open the file to store the results in
  std::ofstream of;
  of.open("mm_timing_res.csv");

  //Initialize cusparse
  
  std::cout << "Created list of matrices will start benchmark: \n";
  //Define the cusparse opaque structures
  cusparseHandle_t handle;
  cusparseStatus_t status = cusparseCreate(&handle);
  
  for(int j = 0; j<files.size();j++)
  {
    try
    {
      run_test_for_matrix(handle,files[j],of);
    }
    catch(int e)
    {
        std::cout << "Caught exception for matrix " << files[j] <<" : " << e << std::endl;
    }
    of.flush();
  }
  of.close();
}

