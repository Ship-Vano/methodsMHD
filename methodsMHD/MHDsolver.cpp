#include"MHDsolver.h"


//const double gam = 5.0 / 3.0; //?
const double gam = 2.0;

//Параметры из вектора состояния
std::vector<double> primitive_vars_from_state(const std::vector<double>& U) {
    std::vector<double> result(8, 0);

    double rho = U[0];
    double mx = U[1];
    double my = U[2];
    double mz = U[3];
    double e = U[4];
    double Bx = U[5];
    double By = U[6];
    double Bz = U[7];

    double vx = mx / rho;
    double vy = my / rho;
    double vz = mz / rho;

    double p = (gam - 1) * (e - 0.5 * rho * (vx * vx + vy * vy + vz * vz)  - 0.5 * (Bx * Bx + By * By + Bz * Bz) );

    result[0] = rho;
    result[1] = vx;
    result[2] = vy;
    result[3] = vz;
    result[4] = p;
    result[5] = Bx;
    result[6] = By;
    result[7] = Bz;

    return result;
}

//Вектор состояния из параметров
std::vector<double> state_from_primitive_vars(double rho, double vx, double vy, double vz, double p, double Bx, double By, double Bz) {
    std::vector<double> U(8,0);

    double mx = rho * vx;
    double my = rho * vy;
    double mz = rho * vz;
    double e = p / (gam - 1) + 0.5 * rho * (vx * vx + vy * vy + vz * vz) + 0.5 * (Bx * Bx + By * By + Bz * Bz);

    U[0] = rho;
    U[1] = mx;
    U[2] = my;
    U[3] = mz;
    U[4] = e;
    U[5] = Bx;
    U[6] = By;
    U[7] = Bz;

    return U;
}

//получить энегию e (то есть 4-ю компоненту состояния U)
double totalenergy(const std::vector<double>& U) {
    return U[4];
}

//определяем скорость быстрой волны
double cfast(const std::vector<double>& U) {
    std::vector<double> components = primitive_vars_from_state(U);
    // 0     1  2   3   4   5  6   7
    //(rho, vx, vy, vz, p, Bx, By, Bz) 
    // *******************************
    //\rho
    double rho = components[0];
    //Bx
    double Bx = components[5];
    //|B|^2
    double BB = Bx * Bx + components[6] * components[6] + components[7] * components[7];
    //|v|^2
    double vv = components[1] * components[1] + components[2] * components[2] + components[3] * components[3];
    //p
    double p = components[4];
    return std::sqrt((gam * p + BB + std::sqrt((gam * p + BB) * (gam * p + BB) - 4 * gam * p * Bx * Bx)) / (2 * rho));
}

//определяем полное давление
double ptotal(const std::vector<double>& U) {
    // 0     1  2   3   4   5  6   7
    //(rho, vx, vy, vz, p, Bx, By, Bz) 
    // *******************************;
    std::vector<double> components = primitive_vars_from_state(U);

    double BB = components[5]* components[5] + components[6] * components[6] + components[7] * components[7];

    return components[4] + BB / 2;
}

//Определяем МГД-поток
std::vector<double> MHD_flux(const std::vector<double>& U) {
    // 0     1  2   3   4   5  6   7
    //(rho, vx, vy, vz, p, Bx, By, Bz) 
    // *******************************;
    std::vector<double> components = primitive_vars_from_state(U);

    double rho = components[0];
    double vx = components[1];
    double vy = components[2];
    double vz = components[3];
    double p = components[4];
    double Bx = components[5];
    double By = components[6];
    double Bz = components[7];

    double BB = Bx * Bx + By * By + Bz * Bz;
    double pT = p + BB / 2;

    std::vector<double> F(8, 0);
    F[0] = rho * vx;
    F[1] = rho * vx * vx + pT - Bx * Bx;
    F[2] = rho * vy * vx - Bx * By;
    F[3] = rho * vz * vx - Bx * Bz;
    F[4] = (totalenergy(U) + pT) * vx - Bx * (vx * Bx + vy * By + vz * Bz);
    F[5] = 0;
    F[6] = By * vx - Bx * vy;
    F[7] = Bz * vx - Bx * vz;

    return F;
}

// Определяем HLL поток F
std::vector<double> HLL_flux(const std::vector<double>& U_L, const std::vector<double>& U_R) {
    // переменные левого и правого края
     // 0    1  2  3  4   5  6   7
    //(rho, u, v, w, p, Bx, By, Bz) 
    // *******************************;
    std::vector<double> componentsL = primitive_vars_from_state(U_L);
    double rhoL = componentsL[0];
    double uL = componentsL[1];
    double vL = componentsL[2];
    double wL = componentsL[3];
    double pL = componentsL[4];
    double BxL = componentsL[5];
    double ByL = componentsL[6];
    double BzL = componentsL[7];
    std::vector<double> componentsR = primitive_vars_from_state(U_R);
    double rhoR = componentsR[0];
    double uR = componentsR[1];
    double vR = componentsR[2];
    double wR = componentsR[3];
    double pR = componentsR[4];
    double BxR = componentsR[5];
    double ByR = componentsR[6];
    double BzR = componentsR[7];

    //быстрые магнитозвуковые скорости на левом и правом концах
    double cfL = cfast(U_L);
    double cfR = cfast(U_R);
    //скорость левого сигнала, рассчитываемая как минимальное значение скорости левого состояния (uL) и быстрой магнитозвуковой скорости (cfL). 
   // double SL = std::min(uL, uR) - std::max(cfL, cfR);
    // Скорость правого сигнала, рассчитываемая как максимальное значение скорости правильного состояния (uR) и быстрой магнитозвуковой скорости (cfR).ы
    double SR = std::max(uL, uR) + std::max(cfL, cfR);
    double SL = -SR;
    // TODO: SL = - SR
    out(U_L);
    std::cout << "SL = " << SL << std::endl;
    if (SL >= 0) {
        return MHD_flux(U_L);
    }
    else if (SR >= 0) {
        return MHD_flux(U_R);
    }
    //(SL <= 0 && SR >= 0)
    else{
        return 1 / (SR - SL) * (SR * MHD_flux(U_L) - SL*MHD_flux(U_R) + SL*SR*(U_R-U_L));
    }

}

//Определяем HLLD поток F*
std::vector<double> HLLD_flux(const std::vector<double>& U_L, const std::vector<double>& U_R) {
    double SL, SR;
    std::vector<double> F_HLLD, F_L, F_R, sF_L, sF_R, ssF_L, ssF_R, sU_r;
    double eL = totalenergy(U_L);
    double eR = totalenergy(U_R);
    // 0    1  2  3  4   5  6   7
    //(rho, u, v, w, p, Bx, By, Bz) 
    // *******************************;
    std::vector<double> componentsL = primitive_vars_from_state(U_L);
    double rhoL = componentsL[0];
    double uL = componentsL[1];
    double vL = componentsL[2];
    double wL = componentsL[3];
    double pL = componentsL[4];
    double BxL = componentsL[5];
    double ByL = componentsL[6];
    double BzL = componentsL[7];
    std::vector<double> componentsR = primitive_vars_from_state(U_R);
    double rhoR = componentsR[0];
    double uR = componentsR[1];
    double vR = componentsR[2];
    double wR = componentsR[3];
    double pR= componentsR[4];
    double BxR = componentsR[5];
    double ByR = componentsR[6];
    double BzR = componentsR[7];

    double Bx = (BxL + BxR) / 2;
    F_L = MHD_flux(U_L);
    F_R = MHD_flux(U_R);
    double cfL = cfast(U_L);
    double cfR = cfast(U_R);
    double pTL = ptotal(U_L);
    double pTR = ptotal(U_R);
    SL = std::min(uL, uR) - std::max(cfL, cfR);
    SR = std::max(uL, uR) + std::max(cfL, cfR);
    double SmuR = SR - uR;
    double SmuL = SL - uL;
    double denSM = SmuR * rhoR - SmuL * rhoL;
    double SM = (SmuR * rhoR * uR - SmuL * rhoL * uL - pTR + pTL) / denSM;
    double suL = SM;
    double ssuL = SM;
    double ssuR = SM;
    double suR = SM;
    double spT = (SmuR * rhoR * pTL - SmuL * rhoL * pTR + rhoL * rhoR * SmuR * SmuL * (uR - uL)) / denSM;
    double spTL = spT;
    double spTR = spT;
    double srhoL = rhoL * SmuL / (SL - SM);
    double srhoR = rhoR * SmuR / (SR - SM);
    double denL = rhoL * SmuL * (SL - SM) - Bx * Bx;
    double denR = rhoR * SmuR * (SR - SM) - Bx * Bx;
    double svL, sByL, swL, sBzL, svR, sByR, swR, sBzR;
    if (std::abs(denL) > 1.0) {
        svL = vL - Bx * ByL * (SM - uL) / denL;
        sByL = ByL * (rhoL * SmuL * SmuL - Bx * Bx) / denL;
        swL = wL - Bx * BzL * (SM - uL) / denL;
        sBzL = BzL * (rhoL * SmuL * SmuL - Bx * Bx) / denL;
    }
    else {
        svL = vL;
        sByL = ByL;
        swL = wL;
        sBzL = BzL;
    }
    if (std::abs(denR) > 1.0) {
        svR = vR - Bx * ByR * (SM - uR) / denR;
        sByR = ByR * (rhoR * SmuR * SmuR - Bx * Bx) / denR;
        swR = wR - Bx * BzR * (SM - uR) / denR;
        sBzR = BzR * (rhoR * SmuR * SmuR - Bx * Bx) / denR;
    }
    else {
        svR = vR;
        sByR = ByR;
        swR = wR;
        sBzR = BzR;
    }
    double vBL = uL * Bx + vL * ByL + wL * BzL;
    double svBL = suL * Bx + svL * sByL + swL * sBzL;
    double seL = (SmuL * eL - pTL * uL + spT * SM + Bx * (vBL - svBL)) / (SL - SM);
    double vBR = uR * Bx + vR * ByR + wR * BzR;
    double svBR = suR * Bx + svR * sByR + swR * sBzR;
    double seR = (SmuR * eR - pTR * uR + spT * SM + Bx * (vBR - svBR)) / (SR - SM);
    double sqrt_srhoL = std::sqrt(srhoL);
    double sqrt_srhoR = std::sqrt(srhoR);
    double sSL = SM - std::abs(Bx) / sqrt_srhoL;
    double sSR = SM + std::abs(Bx) / sqrt_srhoR;
    double ssrhoL = srhoL;
    double ssrhoR = srhoR;
    double sspTL = spTL;
    double sspTR = spTR;
    double sign_Bx = std::copysign(1.0, Bx);
    double sqrt_srhoLsrhoR = sqrt_srhoL * sqrt_srhoR;
    double den = sqrt_srhoL + sqrt_srhoR;
    double ssv = (sqrt_srhoL * svL + sqrt_srhoR * svR + (sByR - sByL) * sign_Bx) / den;
    double ssw = (sqrt_srhoL * swL + sqrt_srhoR * swR + (sBzR - sBzL) * sign_Bx) / den;
    double ssBy = (sqrt_srhoL * sByR + sqrt_srhoR * sByL + sqrt_srhoLsrhoR * (svR - svL) * sign_Bx) / den;
    double ssBz = (sqrt_srhoL * sBzR + sqrt_srhoR * sBzL + sqrt_srhoLsrhoR * (swR - swL) * sign_Bx) / den;
    double ssvL = ssv;
    double ssvR = ssv;
    double sswL = ssw;
    double sswR = ssw;
    double ssByL = ssBy;
    double ssByR = ssBy;
    double ssBzL = ssBz;
    double ssBzR = ssBz;
    double ssvBL = ssuL * Bx + ssvL * ssByL + sswL * ssBzL;
    double ssvBR = ssuR * Bx + ssvR * ssByR + sswR * ssBzR;
    double sseL = seL - sqrt_srhoL * (svBL - ssvBL) * sign_Bx;
    double sseR = seR + sqrt_srhoR * (svBR - ssvBR) * sign_Bx;
    std::vector<double> sU_L{ srhoL, srhoL * suL, srhoL * svL, srhoL * swL, seL, Bx, sByL, sBzL };
    std::vector<double> sU_R{ srhoR, srhoR * suR, srhoR * svR, srhoR * swR, seR, Bx, sByR, sBzR };
    std::vector<double> ssU_L{ ssrhoL, ssrhoL * ssuL, ssrhoL * ssvL, ssrhoL * sswL, sseL, Bx, ssByL, ssBzL };
    std::vector<double> ssU_R{ ssrhoR, ssrhoR * ssuR, ssrhoR * ssvR, ssrhoR * sswR, sseR, Bx, ssByR, ssBzR };
    sF_L = F_L + SL * sU_L - SL * U_L;
    ssF_L = F_L + sSL * ssU_L - (sSL - SL) * sU_L - SL * U_L;
    sF_R = F_R + SR * sU_R - SR * U_R;
    ssF_R = F_R + sSR * ssU_R - (sSR - SR) * sU_R - SL * U_R;
    if (SL > 0) {
        F_HLLD = F_L;
    }
    else if (sSL >= 0) {
        F_HLLD = sF_L;
    }
    else if (SM >= 0) {
        F_HLLD = ssF_L;
    }
    else if (sSR >= 0) {
        F_HLLD = ssF_R;
    }
    else if (SR >= 0) {
        F_HLLD = sF_R;
    }
    else {
        F_HLLD = F_R;
    }
    return F_HLLD;
}

//bool HLLDSolve()
//{
//    string filename = "init1";
//	string path = "OutputData\\" + filename;
//	ofstream fpoints(path);
//	cout << "log[INFO]: Starting HLLD solve..." << endl;
//	cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << endl;
//	//initial values
//    double t_0 = 0.;
//    vector<double> u_0{ 1.,2.,3.,4.,6.,7.,8.};
//    //
//    if (fpoints.is_open())
//	{
//        double t_i = t_0;
//		//vector<double> y_i = u_0;
//		///*fpoints << t_i << endl;
//		//writeVectorToFile(fpoints, y_i);*/
//		//int ind = 0;
//  //      vector<double> y_ipp = u_0;
//		//writeVectorToFile(fpoints, t_i, y_i);
//		//while (abs(T - t_i) >= 1e-8)
//		//{
//		//	y_ipp = y_i + tau * func(t_i, y_i);
//		//	//fpoints << t_i << endl;
//		//	y_i = y_ipp;
//		//	t_i += tau;
//		//	writeVectorToFile(fpoints, t_i, y_i);
//		//}
//		fpoints.close();
//		return true;
//	}
//	else
//		cout << "log[ERROR]: Couldn't open or create a file" << endl;
//	return false;
//}



// #include <iostream>
//#include <vector>
//#include <cmath>
//#include"MHDsolver.h"
//// Constants
//const double gamma = 5.0 / 3.0; // Adiabatic index (specific to MHD)
//const double mu_0 = 1.0; // Magnetic permeability (normalized)
//
//class MHD_solver {
//private:
//	//Grid size and spacing
//	int nx, ny;
//	double dx;
//	//State variables (density, velocity components, pressure, magnetic field components)
//	std::vector<std::vector<double>> state;
//public:
//    //Constructor
//    MHD_solver(int nx, int ny, double dx) : nx(nx), ny(ny), dx(dx) {
//        // Initialize state variables
//        state.resize(static_cast<std::vector<std::vector<double, std::allocator<double>>, std::allocator<std::vector<double, std::allocator<double>>>>::size_type>(nx) * ny, std::vector<double>(8, 0.0)); // Assuming 8 state variables
//        // Initialize other variables and parameters
//        // ...
//    }
//};
////class MHD_Solver {
////public:
////    // Constructor
////    MHD_Solver(int nx, int ny, double dx) : nx(nx), ny(ny), dx(dx) {
////        // Initialize state variables
////        state.resize(static_cast<std::vector<std::vector<double, std::allocator<double>>, std::allocator<std::vector<double, std::allocator<double>>>>::size_type>(nx) * ny, std::vector<double>(8, 0.0)); // Assuming 8 state variables
////        // Initialize other variables and parameters
////        // ...
////    }
////
////    // Method to evolve the system in time
////    void evolve(double dt) {
////        // Loop over each grid point and update the state variables
////        for (int i = 0; i < nx * ny; ++i) {
////            // Compute fluxes using the HLL solver
////            std::vector<double> flux = computeHLLFlux(i);
////
////            // Update state variables using the fluxes and time step
////            updateStateWithFlux(i, flux, dt);
////        }
////    }
////
////private:
////    // Grid size and spacing
////    int nx, ny;
////    double dx;
////
////    // State variables (density, velocity components, pressure, magnetic field components)
////    std::vector<std::vector<double>> state; // 2D array for state variables
////
////    // Method to compute HLL fluxes
////    std::vector<double> computeHLLFlux(int idx) {
////        // Compute left and right states for the Riemann problem
////        std::vector<double> Ul = getStateAtPoint(idx); // Left state
////        std::vector<double> Ur = getStateAtPoint(idx + 1); // Right state
////
////        // Compute and return fluxes using the HLL solver
////        return computeFluxHLL(Ul, Ur);
////    }
////
////    // Method to update state variables with fluxes
////    void updateStateWithFlux(int idx, const std::vector<double>& flux, double dt) {
////        // Update state variables using the fluxes and time step
////        for (int j = 0; j < 8; ++j) { // Assuming 8 state variables
////            state[idx][j] -= dt / dx * (flux[j + 1] - flux[j]);
////        }
////    }
////
////    // Helper method to compute HLL fluxes
////    std::vector<double> computeFluxHLL(const std::vector<double>& Ul, const std::vector<double>& Ur) {
////        // Compute fluxes using the HLL solver
////        // ...
////        std::vector<double> flux;
////        return flux;
////    }
////
////    // Helper method to get state variables at a grid point
////    std::vector<double> getStateAtPoint(int idx) {
////        // Get state variables at grid point idx
////        return state[idx];
////    }
////};
