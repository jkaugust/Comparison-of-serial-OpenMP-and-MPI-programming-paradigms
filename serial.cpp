#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"

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

/*
int buildBins(vector<bin_t>& bins, particle_t* particles, int n)
{
    gridSize = sqrt(n*density);
    binSize = cutoff * 2;  
    binNum = int(gridSize / binSize)+1; // Should be around sqrt(N/2)

    // printf("Grid Size: %.4lf\n",gridSize);
    // printf("Number of Bins: %d*%d\n",binNum,binNum);
    // printf("Bin Size: %.2lf\n",binSize);
    // Increase\Decrease binNum to be something like 2^k?
    
    bins.resize(binNum * binNum);

    for (int i = 0; i < n; i++)
    {
        int x = int(particles[i].x / binSize);
        int y = int(particles[i].y / binSize);
        bins[x*binNum + y].push_back(particles[i]);
    }
    return 0;
}
*/

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
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
    
    //vector<bin_t> particle_bins;
    //vector<bin_t> temp;    
    //buildBins(particle_bins, particles, n);
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    //initialize bins
    for(int i=0; i<totalBins; i++)
    {
      bin_t bin;
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
    for( int step = 0; step < NSTEPS; step++ )
    {
	      navg = 0;
        davg = 0.0;
	      dmin = 1.0;
        //
        //  compute forces
        //
        //clear bins
        for(int i=0;i<totalBins;i++)
        {
          bins[i].particles.clear();
        }
        
        //assign particles to bins
        for(int i=0;i<n;i++)
        {
          int binIndex=find_binIndex(particles[i].x,particles[i].y,binWidth,gridSize);
          particles[i].ax = particles[i].ay = 0;
          bins[binIndex].particles.push_back(i);
        }
        
        //apply force
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
        
        
        /*
        for (int i = 0; i < binNum; i++)
        {
            for (int j = 0; j < binNum; j++)
            {
                bin_t& vec = particle_bins[i*binNum+j];
                for (int k = 0; k < vec.size(); k++)
                    vec[k].ax = vec[k].ay = 0;
                for (int dx = -1; dx <= 1; dx++)   //Search over nearby 8 bins and itself
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        if (i + dx >= 0 && i + dx < binNum && j + dy >= 0 && j + dy < binNum)
                        {
                            bin_t& vec2 = particle_bins[(i+dx) * binNum + j + dy];
                            for (int k = 0; k < vec.size(); k++)
                                for (int l = 0; l < vec2.size(); l++)
                                    apply_force( vec[k], vec2[l], &dmin, &davg, &navg);
                        }
                    }
                }
            }
        }
         */
         
        //
        //  move particles
        //
         for(int i = 0; i < n; i++) 
            move( particles[i] ); 
        /*
        for (int i = 0; i < binNum; i++)
        {
            for(int j = 0; j < binNum; j++)
            {
                bin_t& vec = particle_bins[i * binNum + j];
                int tail = vec.size(), k = 0;
                for(; k < tail; )
                {
                    move( vec[k] );
                    int x = int(vec[k].x / binSize);  //Check the position
                    int y = int(vec[k].y / binSize);
                    if (x == i && y == j)  // Still inside original bin
                        k++;
                    else
                    {
                        temp.push_back(vec[k]);  // Store paricles that have changed bin. 
                        vec[k] = vec[--tail]; //Remove it from the current bin.
                    }
                }
                vec.resize(k);
            }
        }
        
        //
        //  Put the particles in new bin
        //
        for (int i = 0; i < temp.size(); i++)  // Put them into the new bin 
        {
            int x = int(temp[i].x / binSize);
            int y = int(temp[i].y / binSize);
            particle_bins[x*binNum+y].push_back(temp[i]);
        }
        temp.clear();
        */


        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
  }

    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
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
