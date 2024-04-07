#include<iostream>
#include"MHDsolver.h"
#include"FileIO.h"

void brioWu_init(const double& h, std::vector<std::vector<double>>& state)
{
	// x < 0.5
	const double rhoL = 1.;
	const double uL = 0.;
	const double vL = 0.;
	const double wL = 0.;
	const double pL = 1.;
	const double BxL = 0.;
	const double ByL = 1.;
	const double BzL = 0.;

	// x > 0.5
	const double rhoR = 0.125;
	const double uR = 0.;
	const double vR = 0.;
	const double wR = 0.;
	const double pR = 0.1;
	const double BxR = 0.;
	const double ByR = -1.;
	const double BzR = 0.;

	const int num_space_steps = static_cast<int>(1. / h);

	double x_i = 0.;
	for (int i = 0; i <= num_space_steps; ++i)
	{
		std::cout << "x_i = " << x_i << std::endl;
		if (x_i < 0.5) {
			//std::cout << "x_i < 0.5" << std::endl;
			state[i] = state_from_primitive_vars(rhoL, uL, vL, wL, pL, BxL, ByL, BzL);
			out(state[i]);
		}
		else {
			state[i] = state_from_primitive_vars(rhoR, uR, vR, wR, pR, BxR, ByR, BzR);
			out(state[i]);
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
	const double BxL = 0.;
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
	const double BxR = 0.;
	const double ByR = -1.;
	const double BzR = 0.;
	return state_from_primitive_vars(rhoR, uR, vR, wR, pR, BxR, ByR, BzR);
}
int main()
{   
	//Calculation parameters (mesh, steps,....)
	const double x_0 = 0.;
	const double X = 1.;
	const double h = 0.01;
	const double t_0 = 0;
	const double T = 0.1;
	const double tau = 0.001;
	const int num_time_steps = static_cast<int>((T - t_0) / tau);
	const int num_space_steps = static_cast<int>((X-x_0)/ h);
	
	std::vector<std::vector<double>> state_j(num_space_steps+1, std::vector<double>(8,0));
	std::vector<std::vector<double>> state_jpp = state_j;
	std::vector<std::vector<double>> fluxes = state_j;
	//Brio-Wu initializer
	brioWu_init(h, state_j);
	string filename = "BRIOWUtest";
	string path = "OutputData\\" + filename;
	ofstream fpoints(path);
	if (fpoints.is_open())
	{
		double t_i = t_0;
		fpoints << t_0 << std::endl;
		for (int i = 0; i <= num_space_steps; ++i)
		{
			writeVectorToFile(fpoints, state_j[i]);
		}
		fpoints << "--------" << std::endl;
		//System evolution
		for (int j = 0; j <= num_time_steps; ++j) {
			std::cout << "j= " << j << std::endl;
			std::cout << "-----" << std::endl;
			//fluxes[0]= HLL_flux(state_j[0], state_j[0]);
			//fluxes[num_space_steps-1]= HLL_flux(state_j[num_space_steps - 1], state_j[num_space_steps - 1]);
			//for (int i = 1; i < num_space_steps; ++i) {
				//fluxes[i] = HLL_flux(state_j[i-1], state_j[i]);
				//std::cout << i << std::endl;
			//}
			for (int i = 1; i < num_space_steps; ++i) {
				//TODO: внести flux'ы в отдельный вектор и вычислять заранее
				//state_jpp[i] = state_j[i] - tau / h * (fluxes[i] - fluxes[i - 1]);
				state_jpp[i] = state_j[i] - tau / h * (HLL_flux(state_j[i],state_j[i+1])-HLL_flux(state_j[i-1],state_j[i]));
			}
			state_jpp[0] = left_bound();
			state_jpp[num_space_steps] = right_bound();
			state_j = state_jpp;
			t_i += tau;
			fpoints << t_i << std::endl;
			for (int i = 0; i <= num_space_steps; ++i)
			{
				writeVectorToFile(fpoints, state_j[i]);
			}
			fpoints << "--------" << std::endl;
		}
	}
	else
		cout << "log[ERROR]: Couldn't open or create a file" << endl;
	
	/*for(int k = 0; k < num_space_steps; ++k)
		out(state_jpp[k]);*/
	return 0;
}