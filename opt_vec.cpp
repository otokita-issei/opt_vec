#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector3d.h"
#include "opt_vec.h"

int main(void){
    AndoLab::Vector3d <double> r1;
    AndoLab::Vector3d <double> r2;
    AndoLab::Vector3d <double> d;

    r1.set(r_1, th1, phi1, AndoLab::coordinate::Spherical);
    r2.set(r_2, th2, phi2, AndoLab::coordinate::Spherical);
    d = (r2 - r1) / abs(r2 -r1);

    std::string opt_output = "optical_depth/r1_Lat" + std::to_string(Lat_1) + "_Lon" + std::to_string(Lon_1) + "_r2_Lat" + std::to_string(Lat_2) + "_Lon" + std::to_string(Lon_2) + "_Day" + std::to_string(Day) + "_Alt" + std::to_string(z0) + "_" + std::to_string(z_max) + ".dat";
    std::ofstream ofs(opt_output.c_str());

    /*経路1の計算(delta1)*/
    double xi0;
    double t1;
    double delta1;

    xi0 = z0 * KMtoM + R0;      /* xi0 = |r1 + t1d|*/
    t1 = -1.0*r1 % d + std::sqrt((r1 % d)*(r1 % d) + xi0*xi0 - r_1*r_1); /*[m]*/

    delta1 = opt_vertical(z0-1, r1.latitude(), r1.longitude()) * t1 / 1e3;

    ofs << "delta1 = " << delta1 << "\n";

    // double tn_Lm[N+1];
    // double tn_1_Lm[N+1];
    // double tn_Lm_1[N+1];
    // double tn_1_Lm_1[N+1];
    // double Delta [N+1];

    // for(int i = N; i >= 0; i--){
    //     Delta[i] = 0.0;
    // }

    // std::string tnLm_input = "data/Opt_Lat" + std::to_string(Lm) + "_Day" + std::to_string(tn) + ".dat";
    // std::ifstream tnLm(tnLm_input.c_str());
    // std::string tn1_Lm_input = "data/Opt_Lat" + std::to_string(Lm) + "_Day" + std::to_string(tn_1) + ".dat";
    // std::ifstream tn1_Lm(tn1_Lm_input.c_str());
    // std::string tn_Lm1_input = "data/Opt_Lat" + std::to_string(Lm_1) + "_Day" + std::to_string(tn) + ".dat";
    // std::ifstream tn_Lm1(tn_Lm1_input.c_str());
    // std::string tn1_Lm1_input = "data/Opt_Lat" + std::to_string(Lm_1) + "_Day" + std::to_string(tn_1) + ".dat";
    // std::ifstream tn1_Lm1(tn1_Lm1_input.c_str());

    // for(int i = N; i >= 0; i--){
    //     tnLm >> tn_Lm[i];
    //     tn1_Lm >> tn_1_Lm[i];
    //     tn_Lm1 >> tn_Lm_1[i];
    //     tn1_Lm1 >> tn_1_Lm_1[i];
    //     Delta[i] = (tn_1 - t)*(Lm_1 - L)/inc_t/inc_L*tn_Lm[i] + (t - tn)*(Lm_1 - L)/inc_t/inc_L*tn_1_Lm[i] + (tn_1 -t)*(L - Lm)/inc_t/inc_L*tn_Lm_1[i] + (t - tn)*(L - Lm)/inc_t/inc_L*tn_1_Lm_1[i];
    // }

    // delta1 = Delta[z0-1] * t1/1e3;
 
    /*経路2の計算(delta2)*/
    int z_i;
    double xi_i;
    double t_i;
    double delta2 = 0.0;
    for(int i = 0; i < N_z; i++){
        z_i = z0 + i; 
        xi_i = (z_i + 1) * KMtoM + R0;
        t_i = -1.0*r1 % d + std::sqrt((r1 % d)*(r1 % d) + xi_i*xi_i - r_1*r_1);
        delta2 = delta2 + opt_vertical(z_i, (r1 + t_i*d).latitude(), (r1 + t_i*d).longitude()) / (((r1 + t_i*d) % d) / abs(r1 + t_i*d));
        
        ofs << "delta2" << "\n";
        ofs << z_i << "_" << z_i+1 << " " << delta2 << "\n";
    }

    /*経路3の計算(delta3)*/
    double xi_max;
    double t_max;
    double t_R;
    double delta3;

    xi_max = z_max * KMtoM + R0;
    t_max = -1.0*r1 % d + std::sqrt((r1 % d)*(r1 % d) + xi_max*xi_max - r_1*r_1);
    t_R = abs(r2 - r1) - t_max;

    delta3 = opt_vertical(z_max, r2.latitude(), r2.longitude()) * t_R / 1e3;

    ofs << "delta3 = " << delta3 << "\n";

    /*任意の2点の光学的深さ*/
    double Delta;
    Delta = delta1 + delta2 + delta3;

    std::cout << Delta << std::endl;

    ofs << "\n";
    ofs << "Delta = " << Delta << "\n";

    ofs.close();

    return 0;
}