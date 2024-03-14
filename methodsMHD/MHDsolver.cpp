#include <iostream>
#include <vector>
#include <cmath>
#include"MHDsolver.h"
// Constants
const double gamma = 5.0 / 3.0; // Adiabatic index (specific to MHD)
const double mu_0 = 1.0; // Magnetic permeability (normalized)

class MHD_Solver {
public:
    // Constructor
    MHD_Solver(int nx, int ny, double dx) : nx(nx), ny(ny), dx(dx) {
        // Initialize state variables
        state.resize(nx * ny, std::vector<double>(8, 0.0)); // Assuming 8 state variables
        // Initialize other variables and parameters
        // ...
    }

    // Method to evolve the system in time
    void evolve(double dt) {
        // Loop over each grid point and update the state variables
        for (int i = 0; i < nx * ny; ++i) {
            // Compute fluxes using the HLL solver
            std::vector<double> flux = computeHLLFlux(i);

            // Update state variables using the fluxes and time step
            updateStateWithFlux(i, flux, dt);
        }
    }

private:
    // Grid size and spacing
    int nx, ny;
    double dx;

    // State variables (density, velocity components, pressure, magnetic field components)
    std::vector<std::vector<double>> state; // 2D array for state variables

    // Method to compute HLL fluxes
    std::vector<double> computeHLLFlux(int idx) {
        // Compute left and right states for the Riemann problem
        std::vector<double> Ul = getStateAtPoint(idx); // Left state
        std::vector<double> Ur = getStateAtPoint(idx + 1); // Right state

        // Compute and return fluxes using the HLL solver
        return computeFluxHLL(Ul, Ur);
    }

    // Method to update state variables with fluxes
    void updateStateWithFlux(int idx, const std::vector<double>& flux, double dt) {
        // Update state variables using the fluxes and time step
        for (int j = 0; j < 8; ++j) { // Assuming 8 state variables
            state[idx][j] -= dt / dx * (flux[j + 1] - flux[j]);
        }
    }

    // Helper method to compute HLL fluxes
    std::vector<double> computeFluxHLL(const std::vector<double>& Ul, const std::vector<double>& Ur) {
        // Compute fluxes using the HLL solver
        // ...
        std::vector<double> flux;
        return flux;
    }

    // Helper method to get state variables at a grid point
    std::vector<double> getStateAtPoint(int idx) {
        // Get state variables at grid point idx
        return state[idx];
    }
};
