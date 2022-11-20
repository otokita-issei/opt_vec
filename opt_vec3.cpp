#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector3d.h"
#include "opt_vec.h"

int main(void){
    double angle;
    double delta_vertical;
    double delta;
    AndoLab::Vector3d <double> r1;
    AndoLab::Vector3d <double> r2;
    AndoLab::Vector3d <double> r3;
    
    r1.set(r_1, th1, phi1, AndoLab::coordinate::Spherical);
    r2.set(r_2, th2, phi2, AndoLab::coordinate::Spherical);
 
    r3 = r2 - r1;
    r1 = -1.0*r1;
    angle = angle_between(r1, r3);

    std::string opt_output = "optical_depth/r1_Lat" + std::to_string(Lat_1) + "_Lon" + std::to_string(Lon_1) + "_r2_Lat" + std::to_string(Lat_2) + "_Lon" + std::to_string(Lon_2) + "_Day" + std::to_string(Day) + "_Alt" + std::to_string(z0) + "_" + std::to_string(z_max) + ".dat";
    std::ofstream ofs(opt_output.c_str());

    double tn_Lm[N+1];
    double tn_1_Lm[N+1];
    double tn_Lm_1[N+1];
    double tn_1_Lm_1[N+1];
    double Delta[N+1];

        for(int j = N; j >= 0; j--){
            Delta[j] = 0.0;
        }

        std::string tnLm_input = "data/Opt_Lat" + std::to_string(Lm) + "_Day" + std::to_string(tn) + ".dat";
        std::ifstream tnLm(tnLm_input.c_str());

        std::string tn1_Lm_input = "data/Opt_Lat" + std::to_string(Lm) + "_Day" + std::to_string(tn_1) + ".dat";
        std::ifstream tn1_Lm(tn1_Lm_input.c_str());

        std::string tn_Lm1_input = "data/Opt_Lat" + std::to_string(Lm_1) + "_Day" + std::to_string(tn) + ".dat";
        std::ifstream tn_Lm1(tn_Lm1_input.c_str());

        std::string tn1_Lm1_input = "data/Opt_Lat" + std::to_string(Lm_1) + "_Day" + std::to_string(tn_1) + ".dat";
        std::ifstream tn1_Lm1(tn1_Lm1_input.c_str());

        for(int i = N; i >= 0; i--){
            tnLm >> tn_Lm[i];
            tn1_Lm >> tn_1_Lm[i];
            tn_Lm1 >> tn_Lm_1[i];
            tn1_Lm1 >> tn_1_Lm_1[i];
            Delta[i] = (tn_1 - t)*(Lm_1 - L)/inc_t/inc_L*tn_Lm[i] + (t - tn)*(Lm_1 - L)/inc_t/inc_L*tn_1_Lm[i] + (tn_1 -t)*(L - Lm)/inc_t/inc_L*tn_Lm_1[i] + (t - tn)*(L - Lm)/inc_t/inc_L*tn_1_Lm_1[i];
        }

        double z1 = Alt1 * 1.0e-3;
        double z2 = Alt2 * 1.0e-3;

        for(int n = z2; n < z1; n++){
            delta_vertical = delta_vertical + Delta[n];
            ofs << Delta[n]/std::cos(angle) << " " << n << "\n";
        }

    delta = delta_vertical / std::cos(angle);

    // std::cout << "nasukaku = " << angle_between(r1, r3) * Rad2Deg << std::endl;
    // std::cout << "angle = " << angle << std::endl;
    // std::cout << "cos = " << std::cos(angle) << std::endl;
    // std::cout << "delta_vertical = " << delta_vertical << std::endl;  
    std::cout << "delta = " << delta << std::endl;

    ofs << "\n"; 
    ofs << "delta = " << delta << "\n";

    ofs.close();

    return 0; 
   
}