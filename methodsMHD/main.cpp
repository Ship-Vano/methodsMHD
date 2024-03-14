#include<iostream>
#include"MHDsolver.h"
int main()
{   
    // Define grid size and spacing
    int nx = 100, ny = 100;
    double dx = 0.1;

    // Create MHD solver object
    //MHD_Solver solver(nx, ny, dx) {};

    // Time parameters
    double dt = 0.01;
    double total_time = 1.0;

    // Evolve the system in time
    int num_steps = static_cast<int>(total_time / dt);
    for (int step = 0; step < num_steps; ++step) {
   //     solver.evolve(dt);
    }

    // Output results or perform other tasks
	return 0;
}