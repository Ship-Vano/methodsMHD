#include<iostream>
#include"MHDsolver.h"


void brioWu_init(const double& h, std::vector<std::vector<double>>& state)
{
	// x < 0.5
	const double rhoL = 1.;
	const double uL = 0.;
	const double vL = 0.;
	const double wL = 0.;
	const double pL = 1.;
	const double BxL = 0.75;
	const double ByL = 1.;
	const double BzL = 0.;

	// x > 0.5
	const double rhoR = 0.125;
	const double uR = 0.;
	const double vR = 0.;
	const double wR = 0.;
	const double pR = 0.1;
	const double BxR= 0.75;
	const double ByR = -1.;
	const double BzR = 0.;

	const int num_space_steps = static_cast<int>(1 / h);

	double x_i = 0.;
	for (int i = 0; i < num_space_steps; ++i)
	{
		if (x_i < 0.5) {
			state[i] = state_from_primitive_vars(rhoL, uL, vL, wL, pL, BxL, ByL, BzL);
		}
		else {
			state[i] = state_from_primitive_vars(rhoR, uR, vR, wR, pR, BxR, ByR, BzR);
		}
		x_i += h;
	}
}

std::vector<double> left_bound() {
	// x < 0.5
	const double rhoL = 1.;
	const double uL = 0.;
	const double vL = 0.;
	const double wL = 0.;
	const double pL = 1.;
	const double BxL = 0.75;
	const double ByL = 1.;
	const double BzL = 0.;
	return state_from_primitive_vars(rhoL, uL, vL, wL, pL, BxL, ByL, BzL);
}
std::vector<double> right_bound() {
	// x > 0.5
	const double rhoR = 0.125;
	const double uR = 0.;
	const double vR = 0.;
	const double wR = 0.;
	const double pR = 0.1;
	const double BxR = 0.75;
	const double ByR = -1.;
	const double BzR = 0.;
	return state_from_primitive_vars(rhoR, uR, vR, wR, pR, BxR, ByR, BzR);
}
int main()
{   
	//Calculation parameters (mesh, steps,....)
	const double x_0 = 0;
	const double X = 1;
	const double h = 0.001;
	const double t_0 = 0;
	const double T = 0.1;
	const double tau = 0.0001;
	const int num_time_steps = static_cast<int>((T - t_0) / tau);
	const int num_space_steps = static_cast<int>((X-x_0)/ h);
	
	std::vector<std::vector<double>> state_j(num_space_steps, std::vector<double>(8,0));
	std::vector<std::vector<double>> state_jpp = state_j;
	std::vector<std::vector<double>> fluxes = state_j;
	//Brio-Wu initializer
	brioWu_init(h, state_j);

    // Output results or perform other tasks
	return 0;
}