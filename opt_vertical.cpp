#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector3d.h"
#include "opt_vec.h"

double opt_vertical(int z, double Lat, double Lon){
    double lat_L = Lat;      /*任意の緯度*/
    int L_m = lat_L/10; 
    int Lm = L_m*10;  /*切り捨てたlat_L*/
    int Lm_1 = Lm +10;
    double L = lat_L;
    int N = z_km;

    double tn_Lm[N+1];
    double tn_1_Lm[N+1];
    double tn_Lm_1[N+1];
    double tn_1_Lm_1[N+1];
    double Delta[N+1];

    for(int i = N; i >= 0; i--){
        Delta[i] = 0.0;
    }

    std::string tnLm_input = "data/Opt_Lat" + std::to_string(Lm) + "_Day" + std::to_string(tn) + ".dat";
    std::ifstream tnLm(tnLm_input.c_str());    std::string tn1_Lm_input = "data/Opt_Lat" + std::to_string(Lm) + "_Day" + std::to_string(tn_1) + ".dat";
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

        return Delta[z];
}