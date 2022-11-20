/*opt_vec2.cpp*/

/*変化させるパラメタ*/
/*1つめの任意の点*/
constexpr double Lat1 {42};
constexpr double Lon1 {0};
constexpr double Alt1 {54.9e3}; /*高度[m]*/

/*2つめの任意の点*/
constexpr double Lat2 {85};
constexpr double Lon2 {0};
constexpr double Alt2 {60.1e3}; /*Alt1 < Alt2*/

constexpr double day_t {172};     /*任意の日にち*/

/*固定するパラメタ*/

constexpr double LONGITUDE_HIMAWARI {140.7}; /*気象衛生ひまわりの経度*/
constexpr double Deg2Rad {M_PI/180.0};  /*DegをRadに変換*/
constexpr double Rad2Deg {180.0/M_PI};  /*RadをDegに変換*/
constexpr double MtoKM {1.0e-3};    /*mをkmに変換*/
constexpr double KMtoM {1.0e3};     /*kmをmに変換*/
constexpr double R0 {6370.0e3};     /*地球の半径[m]*/

/*1つめの任意の点の球座標*/
constexpr double r_1 {R0 + Alt1};
constexpr double th1 {(90 - Lat1) * Deg2Rad};
constexpr double phi1 {(Lon1 - LONGITUDE_HIMAWARI) * Deg2Rad};

/*2つめの任意の点の球座標*/
constexpr double r_2 {R0 + Alt2};
constexpr double th2 {(90 - Lat2) * Deg2Rad};
constexpr double phi2 {(Lon2 - LONGITUDE_HIMAWARI) * Deg2Rad};

constexpr int z0 = Alt1*MtoKM + 1;
constexpr int z_max = Alt2*MtoKM;
constexpr int N_z = z_max - z0 + 1;

constexpr int Lat_1 = Lat1; /*opt_output用*/
constexpr int Lat_2 = Lat2;
constexpr int Lon_1 = Lon1;
constexpr int Lon_2 = Lon2;
constexpr int Day = day_t;

/*opt_vertical*/
double opt_vertical(int z, double Lat, double Lon);
constexpr double lat_L = Lat1;      /*任意の緯度*/
constexpr double z_km {100};    /*高度*/
constexpr double t0 = 100;      /*基準の日にち*/
constexpr double L0 = 40;       /*基準の緯度*/
constexpr double inc_t = 10;    /*日にちの増加量*/
constexpr double inc_L = 10;    /*緯度の増加量*/

constexpr int t_n = day_t/10; 
constexpr int tn = t_n*10;  /*切り捨てたday_t*/
constexpr int tn_1 = tn + 10;

constexpr int L_m = lat_L/10; 
constexpr int Lm = L_m*10;  /*切り捨てたlat_L*/
constexpr int Lm_1 = Lm +10;

constexpr double t = day_t;
constexpr double L = lat_L;
//constexpr double h_max = z_km;
constexpr int N = z_km;

// double tn_Lm[N+1];
// double tn_1_Lm[N+1];
// double tn_Lm_1[N+1];
// double tn_1_Lm_1[N+1];
// double Delta[N+1];

// std::string tnLm_input = "data/Opt_Lat" + std::to_string(Lm) + "_Day" + std::to_string(tn) + ".dat";
// std::ifstream tnLm(tnLm_input.c_str());

// std::string tn1_Lm_input = "data/Opt_Lat" + std::to_string(Lm) + "_Day" + std::to_string(tn_1) + ".dat";
// std::ifstream tn1_Lm(tn1_Lm_input.c_str());

// std::string tn_Lm1_input = "data/Opt_Lat" + std::to_string(Lm_1) + "_Day" + std::to_string(tn) + ".dat";
// std::ifstream tn_Lm1(tn_Lm1_input.c_str());

// std::string tn1_Lm1_input = "data/Opt_Lat" + std::to_string(Lm_1) + "_Day" + std::to_string(tn_1) + ".dat";
// std::ifstream tn1_Lm1(tn1_Lm1_input.c_str());

// std::ifstream tnLm("Opt_Lat40_Day150.dat");
// std::ifstream tn1_Lm("Opt_Lat40_Day160.dat");
// std::ifstream tn_Lm1("Opt_Lat50_Day150.dat");
// std::ifstream tn1_Lm1("Opt_Lat50_Day160.dat");