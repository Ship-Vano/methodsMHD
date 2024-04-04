#pragma once
#include<cmath>
#include<vector>
#include"LinOp.h"
#include"FileIO.h";

std::vector<double> primitive_vars_from_state(const std::vector<double>& U);

std::vector<double> state_from_primitive_vars(double rho, double vx, double vy, double vz, double p, double Bx, double By, double Bz);

double cfast(const std::vector<double>& U);

std::vector<double> MHD_flux(const std::vector<double>& U);

std::vector<double> HLL_flux(const std::vector<double>& U_L, const std::vector<double>& U_R);
