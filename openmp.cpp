#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
#include "omp.h"
using std::vector;

#define cutoff 0.01    
#define density 0.0005

//Find Bins for storing the particles which are relatively close to each other
int find_binIndex(double x_direction, double y_direction, double bin_size, double grid_block_size)
{
  int bin_x,bin_y;
  //find x coordinate of bin

  bin_x=int(x_direction/bin_size);
  
  //find y coordinate of bin
  bin_y=int(y_direction/bin_size);
  
  return(bin_y*(grid_block_size/bin_size)+bin_x);
}


//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    double binWidth,gridSize;
    int binNum;
    gridSize = sqrt(n*density);
    binWidth = sqrt(3*density);  
    binNum = int(gridSize / binWidth)+1;
    //binNum = int(gridSize / binSize);
    binWidth = gridSize/binNum;
    
    int totalBins=binNum*binNum;
    bin_t bins[totalBins];
    omp_lock_t locks[totalBins];
    //vector<bin_t> particle_bins;
    //vector<bin_t> temp;
    //buildBins(particle_bins, particles, n);
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    //initialize bins
    #pragma omp for
    for(int i=0; i<totalBins; i++)
    {
      bin_t bin;
      omp_init_lock(locks+i);
      int x_bin=i%binNum;
      int y_bin=int(i/binNum);
      if(x_bin==0||x_bin==(binNum-1)||y_bin==0||y_bin==(binNum-1))
      {
      //check neighbours
        for(int x=-1;x<=1;x++)
        {
          for(int y=-1;y<=1;y++)
          {
            if(x_bin+x>=0 && y_bin+y>=0 && x_bin+x<binNum && y_bin+y<binNum)
            {
              bin.neighbours.push_back(i+x+binNum*y);
            }
          }
        }
      }else
      {
        for(int x=-1;x<=1;x++)
        {
          for(int y=-1;y<=1;y++)
          {
            bin.neighbours.push_back(i+x+binNum*y);
          }
        }
      }
      bins[i]=bin;
    }
    
    //Simulate N time steps
    #pragma omp parallel private(dmin)
    {
      numthreads=omp_get_num_threads();
    for( int step = 0; step < NSTEPS; step++ )
    {
	      navg = 0;
        davg = 0.0;
	      dmin = 1.0;
        //
        //  compute forces
        //
        //clear bins
        #pragma omp for
        for(int i=0;i<totalBins;i++)
        {
          bins[i].particles.clear();
        }
        //assign particles to bins
        #pragma omp for
        for(int i=0;i<n;i++)
        {
          int binIndex=find_binIndex(particles[i].x,particles[i].y,binWidth,gridSize);
          particles[i].ax = particles[i].ay = 0;
          omp_set_lock(locks+binIndex);
          bins[binIndex].particles.push_back(i);
          omp_unset_lock(locks+binIndex);
        }
        //apply force
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for(int i=0;i<totalBins;i++)
        {
          for(std::vector<int>::size_type j=0;j!=bins[i].neighbours.size();j++)
          {
            for(std::vector<int>::size_type k=0;k!=bins[i].particles.size();k++)
            {
              for(std::vector<int>::size_type l=0;l!=bins[bins[i].neighbours[j]].particles.size();l++)
              {
                apply_force(particles[k],particles[l],&dmin,&davg,&navg);
              }
            }
          }
        }
        
        //
        //  move particles
        //
        #pragma omp for
         for(int i = 0; i < n; i++) 
            move( particles[i] );
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
        #pragma omp barrier
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
