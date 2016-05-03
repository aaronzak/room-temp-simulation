//midterm.c
// https://people.sc.fsu.edu/~jburkardt/c_src/heated_plate/heated_plate.html
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
//#include <omp.h>


// Prototypes of Helper functions:
double x(int i); 
double y(int i);                   // x-coordinate in the [a, b] interval corresponding to index i.
double exact_solution(double x, double y);    // Definition of the exact solution.
double S(double x, double y);                 // Definition of the source term.

// Givens of the problem:
const double a=0.0; const double b=1.0; // Interval [a, b] chosen.
const double c=0.0; const double d=2.0;
const int m=40;                         // Number of discretization points chosen.
const int n=80;
// Step sizes:
const double dx=(b-a)/(m-1);
const double dy=(d-c)/(m-1);


// We need 2 arrays, Un and Unp1, for computing the solution:
double Un[m][n];
double Unp1[m][n];

int main()
{
    int thread_count;
    std::cin>>thread_count;


    // index:
    int i;
    int j;

    double time1;
    double time2;
    double timeD;
    // Variables used to establish the convergence of the Jacobi iterations:
    double iteration_error=1.0;
    double tolerance=1E-15;
    int Max_Iter=1000000;

    // Initialize Un arbitrarily to 0:
    #pragma omp parallel for num_threads(thread_count) shared ( Un ) private ( i, j )
    for(i=1; i<m-1; i++) {
        for(j=1;j<n-1;j++){
            Un[i][j]=0.0;
        }}

    // Impose the left and right boundary conditions:

    Un[0][0]   =  exact_solution(a,c);
    Un[m-1][0] =  exact_solution(b,c);
    Un[0][n-1] =  exact_solution(a,d);
    Un[m-1][n-1]=  exact_solution(b,d);

    // Initialize the interation counter:
    int iteration_count=0;
    
         
        

    // Iterate until steady-state is reached.
    // Otherwise stops at Max_Iter iterations (to avoid infinite loops).
    //
    // Recall that we say that the steady-state is reached when the maximum difference between
    // two iterates is less than or equal to the tolerance, i.e. max|Unp1-Un| <= tolerance.

    time1=omp_get_wtime();

    while(iteration_error > tolerance && iteration_count < Max_Iter){
        iteration_count++; // if(iteration_count % 1000 == 0) std::cout<<"iteration " << iteration_count << std::endl;

        // Treat the left and right boundary conditions:

        Unp1[0][0]   = exact_solution(a,c); 
        Unp1[m-1][0] = exact_solution(b,c);
         Unp1[0][n-1] = exact_solution(a,d);
     Unp1[m-1][n-1]= exact_solution(b,d);

        // Treat the interior points using the update formula:

       #pragma omp parallel for num_threads(thread_count) shared ( Un, Unp1 ) private ( i, j )
        for(i=1; i< m-1; i++){
            for(j=1; j< n-1;j++){
                #pragma omp critical
            Unp1[i][j]=.25*( Un[i+1][j] + Un[i-1][j] + Un[i][j-1] + Un[i][j+1]-  dx*dy*S(x(i),y(j)) );
        }}
        

        // Compute the maximum error between 2 iterates to establish whether or not
        // steady-state is reached:
          iteration_error=0.0;
          #pragma omp parallel for num_threads(thread_count) shared ( Un, Unp1, iteration_error ) private ( i, j )
        for(i=0; i< m; i++){
            for(j=1;j<n-1;j++){
            double local_iteration_error = fabs(Unp1[i][j] - Un[i][j]);
            if (local_iteration_error > iteration_error) {
                #pragma omp critical
                iteration_error = local_iteration_error;
            }
        }
        }
        

        // Prepare for the next iteration:
   #pragma omp parallel for num_threads(thread_count) shared ( Un, Unp1 ) private ( i, j )
            for(i=0; i< m; i++){
            for(j=1;j<n-1;j++){
            Un[i][j]=Unp1[i][j];
        }
        }
        

//        if(iteration_count % 1000 == 0) std::cout<< "The error between two iterates is " << iteration_error << std::endl;
    }

    // Compute the maximum error between the computed and exact solutions:
        double solution_error=0.0;
        #pragma omp parallel for num_threads(thread_count) shared ( Un,Unp1, solution_error ) private ( i, j )
    for(i=0; i< m; i++){
        for(j=0; j<m;j++){
            double unp1var= Unp1[i][j];
        double local_solution_error=fabs(unp1var - exact_solution(x(i),y(j)) );
        if (local_solution_error > solution_error) {
             #pragma omp critical
            solution_error = local_solution_error;}
    }
    }

    time2=omp_get_wtime();
    timeD=time2-time1;

for(i=0; i< m; i++){
            for(j=1;j<n-1;j++){
           std::cout<<Un[i][j]<<std::endl;
        }
        }
    // Output:
    printf("\n\n");
    std::cout<< "-------------------------------------------------------"               << std::endl;
    std::cout<< "SUMMARY:"                                                 << std::endl << std::endl;
    std::cout<< "The error between two iterates is "    << iteration_error << std::endl << std::endl;
    std::cout<< "The maximum error in the solution is " << solution_error               << std::endl;
    std::cout<< "iteration count" <<iteration_count<<std::endl;  
    std::cout<< "time taken " <<timeD<<std::endl;

    std::cout<< "-------------------------------------------------------"  << std::endl << std::endl;

 
    return 0;
}

// Helper functions:

double x(int i){
    return a+i*dx;
}
double y(int i){
    return c+i*dy;
}

double exact_solution(double x, double y){
    return sin(2*M_PI*x)+cos(2*M_PI*y);
}

double S(double x,double y)
{
    return -8*M_PI*M_PI*exact_solution(x,y);
}