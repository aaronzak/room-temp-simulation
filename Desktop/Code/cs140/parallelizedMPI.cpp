//midtermcoutn-.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <iostream>


// Prototypes of Helper functions:
double x(int i); 
double y(int i);                   // x-coordinate in the [a, b] interval corresponding to index i.
double exact_solution(double x, double y);    // Definition of the exact solution.
double S(double x, double y);                 // Definition of the source term.

// Givens of the problem:
const double a=0.0; const double b=1.0; // Interval [a, b] chosen.
const double c=0.0; const double d=1.0;
const int m=40;                         // Number of discretization points chosen.
const int n=40;
// Step sizes:
const double dx=(b-a)/(m-1);
const double dy=(d-c)/(n-1);
const double dx2=dx*dx;
const double dy2=dy*dy;
const double dxdyd= 1/(dx*dx+dy*dy);
// We need 2 arrays, Un and Unp1, for computing the solution:
double Un[m][n];
double Unp1[m][n];

int main()
{


  int rank, processes, process;
   MPI_Status status;
   int master=0;
   int tag=123;

   int thread_count=1;
    MPI_Init(NULL,NULL);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // index:
    int i;
    int j;
    MPI_Request requests[processes];
    double time1;
    double time2;
    double timeD;
    // Variables used to establish the convergence of the Jacobi iterations:
    double iteration_error=1.0;
    double tolerance=1E-15;
    int Max_Iter=1000000/processes;
    if(rank==0){
    double d_sqrt = sqrt(processes);
    int i_sqrt = d_sqrt;
    if ( d_sqrt != i_sqrt){

      fprintf(stderr,"Number of processes must be perfect square\n");
      MPI_Finalize();
    return 0;
    }
  }

    // Initialize Un arbitrarily to 0:
    //     #pragma omp parallel for num_threads(thread_count) shared ( Un ) private ( i, j )
  int sqpr=sqrt(processes);
  int m_iter=m/sqpr;
  int n_iter=n/sqpr; 
// create a 2d array for each processes and set the arrays to 0
  double local_Un[m_iter][n_iter];
  double local_Unp1[m_iter][n_iter];


//MPI: seperate the work between threads into processes number of 2d arrays
  //  0  1   2   3
  //  4  5   6   7
  //  8  9   10  11
  //  12 13  14  15


  

    for(i=0; i<m_iter; i++) {
        for(j=0; j<n_iter; j++){
    local_Un[i][j]=0;
        }}
    
          
    // Impose the left and right boundary conditions:

        int local_column = rank%sqpr;
        int local_row = rank/sqpr;

    if(local_row==0){
       for(i=0; i<n_iter; i++) {
        local_Un[0][i]= exact_solution(x(i+(local_column)*n_iter),c);
       
       }
    }

    if(local_row==(sqpr-1)){
      for(i=0; i<n_iter; i++) {
        local_Un[m_iter-1][i]= exact_solution(x(i+(local_column)*n_iter),d);
       }
    }

    if(local_column==0){
      for(j=0;j<m_iter;j++){
        local_Un[j][0]= exact_solution(a,y(j+(local_row)*m_iter));
      }
    }
    if(local_column==(sqpr-1)) {
      for(j=0;j<m_iter;j++){
        local_Un[j][n_iter-1]= exact_solution(b,y(j+(local_row)*m_iter));
      }
    }


    /*

    for(i=0;i<m_iter;i++){
      local_Un[i][0]= exact_solution(x(0+(local_column)*m_iter),y(i+(local_row)*m_iter));
      local_Un[i][n_iter-1]= exact_solution(x(n_iter-1+(local_column)*m_iter),y(i+(local_row)*m_iter));


  }

      for(j=0;j<m_iter;j++){
        local_Un[0][j]= exact_solution(x(j+(local_column)*m_iter),y(0+(local_row)*m_iter));
        local_Un[n_iter-1][j]= exact_solution(x(j+(local_column)*m_iter),y(n_iter-1+(local_row)*m_iter));

      }
    

    /*for(i=0;i<m_iter;i++){
      for(j=0;j<n_iter;j++){
        printf("[%d][%d]= %f",i+(rank/sqpr)*m_iter,j+(rank%sqpr)*m_iter,local_Un[i][j]);
      }
      printf("\n\n");
    }*/

    // Initialize the interation counter:
    int iteration_count;
    iteration_count=0;
    
  
    // Iterate until steady-state is reached.
    // Otherwise stops at Max_Iter iterations (to avoid infinite loops).
    //
    // Recall that we say that the steady-state is reached when the maximum difference between
    // two iterates is less than or equal to the tolerance, i.e. max|Unp1-Un| <= tolerance.
    

       time1=MPI_Wtime();




  MPI_Barrier(MPI_COMM_WORLD);



       
       // {
    while(iteration_error> tolerance && iteration_count < Max_Iter){
      {      
       //iteration_count++;
	/*	if(rank==0){
	  iteration_count++;
	  if(iteration_count>10){
	    MPI_Bcast(&iteration_count,1,MPI_INT,0,MPI_COMM_WORLD);
	  }
	  }*/
	iteration_count++;


    double upToDown[n_iter];
    double downToUp[n_iter];
    double leftToRight[m_iter];
    double rightToLeft[m_iter];
       // if(iteration_count % 1000 == 0) std::cout<<"iteration " << iteration_count << std::endl;
//broadcast iteration count to all processes
       

       // Treat the left and right boundary conditions:

       if(local_row==0){
       for(i=0; i<n_iter; i++) {
        local_Unp1[0][i]= exact_solution(x(i+(local_column)*n_iter),c);
       }

       }

    if(local_row==(sqpr-1)){
      for(i=0; i<n_iter; i++) {
        local_Unp1[m_iter-1][i]= exact_solution(x(i+(local_column)*n_iter),d);
       }
    }
    if((local_column)!=(sqpr-1)){
      for(i=0;i<n_iter;i++){
        leftToRight[i]=  exact_solution(x(m_iter-1+local_column*m_iter),y(i+local_row*m_iter));
      }  
     }
    if(local_column==0){
      for(j=0;j<m_iter;j++){
        local_Unp1[j][0]= exact_solution(a,y(j+(local_row)*m_iter));
      }
    }
    if((local_row)!=(sqpr-1)){
        for(i=0;i<n_iter;i++){
          upToDown[i]= exact_solution(x(i+local_column*m_iter),y(m_iter-1+local_row*m_iter));
       }
     }
    if(local_column==(sqpr-1)) {
      for(j=0;j<m_iter;j++){
        local_Unp1[j][n_iter-1]= exact_solution(b,y(j+(local_row)*n_iter));
      }
    }
     


    for(i=0;i<m_iter;i++){
      local_Unp1[i][0]= exact_solution(x(0+(local_column)*m_iter),y(i+(local_row)*m_iter));
      local_Unp1[i][n_iter-1]= exact_solution(x(n_iter-1+(local_column)*m_iter),y(i+(local_row)*m_iter));


  }

      for(j=0;j<m_iter;j++){
        local_Unp1[0][j]= exact_solution(x(j+(local_column)*m_iter),y(0+(local_row)*m_iter));
        local_Unp1[n_iter-1][j]= exact_solution(x(j+(local_column)*m_iter),y(n_iter-1+(local_row)*m_iter));

      }

    

   /* printf("rank is %d",rank);
    for(i=0;i<m_iter;i++){
      for(j=0;j<n_iter;j++){
        printf("[%d][%d]= %f",i+(rank/sqpr)*m_iter,j+(rank%sqpr)*m_iter,local_Un[i][j]);
      }
      printf("\n\n");
    }
*/


  
   



        // Treat the interior points using the update formula:
	
//parallelize the update for loop and use scatter to spread the solutions
   

   
       

     for(i=1; i< m_iter-1; i++){
	  
            for(j=1; j< n_iter-1;j++){
	      local_Unp1[j][i]= 0.5*(dy2*(local_Un[j+1][i]+local_Un[j-1][i])+dx2*(local_Un[j][i+1]+local_Un[j][i-1]) - S(x(i+(local_column)*m_iter),y(j+(local_row)*n_iter))*dx2*dy2)*dxdyd;
  }
}

     //send from up to Down
  //all 2d besides the top row and grap the ghost from above
/*

 if((local_row)!=(sqpr-1)){
      for(i=0;i<n_iter;i++){
       // upToDown[i]= local_Unp1[m_iter-1][i];
      }

       MPI_Send(upToDown,n_iter,MPI_DOUBLE,rank+sqpr,100,MPI_COMM_WORLD);  
     }


     if((local_row)!=0){
       MPI_Recv(upToDown,n_iter,MPI_DOUBLE,rank-sqpr,100,MPI_COMM_WORLD,&status);
       
       for(i=1;i<n_iter-1;i++){
	 local_Unp1[0][i]= 0.5*(dy2*(local_Un[0][i+1]+local_Un[0][i-1])+dx2*(local_Un[1][i]+upToDown[i]) - S(x(i+(local_column)*n_iter),y(0+(local_row)*m_iter))*dx2*dy2)*dxdyd;
	 downToUp[i]= local_Unp1[0][i];
       }
       downToUp[0]=local_Unp1[0][0];
       downToUp[n_iter-1]=local_Unp1[0][n_iter-1];
       MPI_Send(downToUp,n_iter,MPI_DOUBLE,rank-sqpr,100,MPI_COMM_WORLD);
     
	 
   
	 
     }


     //all arrays besides bottom get the ghost array from their underprocess
     if((local_row)!=(sqpr-1)){
      MPI_Recv(downToUp,n_iter,MPI_DOUBLE,rank+sqpr,100,MPI_COMM_WORLD,&status);
      
      
        for(i=1;i<n_iter-1;i++){
   local_Unp1[m_iter-1][i]= 0.5*(dy2*(local_Un[m_iter-1][i+1]+local_Un[m_iter-1][i-1])+dx2*(local_Un[m_iter-2][i]+downToUp[i]) - S(x(i+(local_column)*n_iter),y(m_iter-1+(local_row)*m_iter))*dx2*dy2)*dxdyd;
       }

     
     }




////////////left to right

      
    


     if((local_column)!=(sqpr-1)){

      MPI_Send(leftToRight,n_iter,MPI_DOUBLE,rank+1,300,MPI_COMM_WORLD);

     }

//all 2d besides the left side so we grab the ghopst from the left
     if((local_column)!=0){

      MPI_Recv(leftToRight,n_iter,MPI_DOUBLE,rank-1,300,MPI_COMM_WORLD,&status);
        for(i=1;i<n_iter-1;i++){
   local_Unp1[i][0]= 0.5*(dy2*(local_Un[i][1]+leftToRight[i])+dx2*(local_Un[i+1][0]+local_Un[i-1][0]) - S(x(0+(local_column)*n_iter),y(i+(local_row)*m_iter))*dx2*dy2)*dxdyd;
   rightToLeft[i]=local_Unp1[i][0];
  //    rightToLeft[i]=exact_solution(x(0+local_column*m_iter),y(i+local_row*m_iter));

       }
       rightToLeft[0]=local_Unp1[0][0];
       rightToLeft[n_iter-1]=local_Unp1[n_iter-1][0];
       MPI_Send(rightToLeft,n_iter,MPI_DOUBLE,rank-1,400,MPI_COMM_WORLD);

       

     }

     //all 2d array besides the right
     if((local_column)!=(sqpr-1)){

      MPI_Recv(rightToLeft,n_iter,MPI_DOUBLE,rank+1,400,MPI_COMM_WORLD,&status);
      for(i=1;i<n_iter-1;i++){
        local_Unp1[i][m_iter-1]= 0.5*(dy2*(local_Un[i][m_iter-2]+rightToLeft[i])+dx2*(local_Un[i+1][m_iter-1]+local_Un[i-1][m_iter-1]) - S(x(m_iter-1+(local_column)*n_iter),y(i+(local_row)*m_iter))*dx2*dy2)*dxdyd;
      }
           
     }

     ///////////




     //work on the bot lef corner
     if(local_column!=0&& local_row!=(sqpr-1)){
        local_Unp1[n_iter-1][0]= 0.5*(dy2*(local_Un[n_iter-1][1]+leftToRight[n_iter-1])+dx2*(local_Un[n_iter-2][0]+downToUp[0]) - S(x(0+(local_column)*n_iter),y(n_iter-1+(local_row)*m_iter))*dx2*dy2)*dxdyd;
        local_Unp1[n_iter-1][0]=exact_solution(x(0+(local_column)*m_iter),y(n_iter-1+(local_row)*m_iter));
     }

     //bottomRight
     if(local_column!=(sqpr-1) && local_row!=(sqpr-1)){
      local_Unp1[n_iter-1][m_iter-1]= 0.5*(dy2*(local_Un[n_iter-1][m_iter-2]+rightToLeft[n_iter-1])+dx2*(local_Un[n_iter-2][m_iter-1]+downToUp[m_iter-1]) - S(x(m_iter-1+(local_column)*n_iter),y(n_iter-1+(local_row)*m_iter))*dx2*dy2)*dxdyd;
      local_Unp1[n_iter-1][n_iter-1] = exact_solution(x(n_iter-1+(local_column)*m_iter),y(n_iter-1+(local_row)*m_iter));
     }

     //top left
     if(local_column!=0 && local_row!=0){
        local_Unp1[0][0]= 0.5*(dy2*(local_Un[0][1]+leftToRight[0])+dx2*(local_Un[1][0]+upToDown[0]) - S(x(0+(local_column)*n_iter),y(0+(local_row)*m_iter))*dx2*dy2)*dxdyd;
        local_Unp1[0][0] = exact_solution(x(0+(local_column)*m_iter),y(0+(local_row)*m_iter));

     }

     //top right
     if(local_column!=(sqpr-1)&&local_row!=0){
        local_Unp1[0][m_iter-1]= 0.5*(dy2*(local_Un[0][m_iter-2]+rightToLeft[0])+dx2*(local_Un[1][m_iter-1]+downToUp[m_iter-1]) - S(x(m_iter-1+(local_column)*n_iter),y(0+(local_row)*m_iter))*dx2*dy2)*dxdyd;
        local_Unp1[0][n_iter-1] = exact_solution(x(0+(local_column)*m_iter),y(n_iter-1+(local_row)*m_iter));
     }

//MPI_Waitall(processes,requests,MPI_STATUSES_IGNORE);
      
	//	
        // Compute the maximum error between 2 iterates to establish whether or not
        // steady-state is reached:

     

*/


	double my_iteration_error=0.0;
	  for(i=0; i< m_iter-1; i++){
            for(j=0;j<n_iter-1;j++){
            double local_iteration_error = fabs(local_Unp1[i][j] - local_Un[i][j]);
	  
	    if (local_iteration_error > my_iteration_error) {
	      // #pragma omp critical
        my_iteration_error = local_iteration_error;

        	       
	    } 
	    }
	  }



    if(rank==0){
      
      
      if(iteration_error > my_iteration_error){
          iteration_error=my_iteration_error;
        }
      for(i=1;i<processes;i++){
        
        MPI_Recv(&my_iteration_error,1,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);

        if(iteration_error > my_iteration_error){
          iteration_error=my_iteration_error;
        }
      }
    }
    else{
    MPI_Send(&my_iteration_error,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
        }

	

    MPI_Bcast(&iteration_error,1,MPI_DOUBLE,0,MPI_COMM_WORLD);


            MPI_Barrier(MPI_COMM_WORLD);
    
    //Do I just make everything Isend then call wait??


        // Prepare for the next iteration:
        //Wait for all the updates to finish then we can start preperation
            for(i=0; i< m_iter; i++){
            for(j=0;j<n_iter;j++){
            local_Un[i][j]=local_Unp1[i][j];
	    }
	    }
	/*	*/   
	
//        if(iteration_count % 1000 == 0) std::cout<< "The error between two iterates is " << iteration_error << std::endl;
      }
      // std::cout<<Unp1[5][5]<<","<<exact_solution(x(5),y(5))<<std::endl;
    } //Done with While loop


                //MPI_Barrier(MPI_COMM_WORLD);

    time2=MPI_Wtime();
    






    int ii;
    int jj;

    int counter=0;

// Compute the maximum error between the computed and exact solutions:
        double solution_error=0.0;
        double my_solution_error=0.0;

    for(i=1; i< m_iter-1; i++){
        for(j=1; j<n_iter-1;j++){
	 
	 
	
	  double local_solution_error=fabs(local_Unp1[j][i] - exact_solution(x(i+(local_column)*m_iter),y(j+(local_row)*n_iter)));
       if(rank==10){  std::cout<<"error of rank "<<rank<<" is "<<local_solution_error<<" "<<j<<i<<std::endl;
          }
	if (local_solution_error > my_solution_error) {
//	
            my_solution_error = local_solution_error;}
	
	}
    }



if(rank==0){
  std::cout<<"checking solution error in master "<<my_solution_error<<std::endl;
      solution_error=my_solution_error;
      for(i=1;i<processes;i++){
        MPI_Recv(&my_solution_error,1,MPI_DOUBLE,i,2,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
	//std::cout<<"Rank "<<rank<<"error"<<my_solution_error<<std::endl;
	if(my_solution_error > solution_error){

    

          solution_error=my_solution_error;
        }
      }
    }
    else{
    MPI_Send(&my_solution_error,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
        }


	
    
    


    // time2=omp_get_wtime();
    timeD=time2-time1;

 double diff;

 if(rank==0){
    // Output:
    printf("\n\n");
    std::cout<< "-------------------------------------------------------"               << std::endl;
    std::cout<< "SUMMARY:"                                                 << std::endl << std::endl;
    std::cout<< "The error between two iterates is "    << iteration_error << std::endl << std::endl;
    std::cout<< "The maximum error in the solution is " << solution_error <<std::endl;
    std::cout<< "iteration count" <<iteration_count<<" per process"<<std::endl;  
    std::cout<< "time taken " <<timeD<<std::endl;
    std::cout<< "-------------------------------------------------------"  << std::endl << std::endl;
 }
/*
for(i=0; i< m; i++){
  for(j=1;j<n-1;j++){
    std::cout<<counter<<","<<Unp1[i][j]<<std::endl;
    counter++;
  }
  } */
 MPI_Finalize();
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
    return sin(2*M_PI*x)*cos(2*M_PI*y);
}

double S(double x,double y)
{
  return -8*M_PI*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y);
}
