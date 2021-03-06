const double rho0 = 2300.0;
const double p0 = 1.0e5;
const double c0 = 1500.0;
const double g = 10.0;
const double Kh = 4.33e-6;
const double kB = 1.38e-23;
const double Rg = 8.314;
const double T = 1000;
const double E0 = 5.1e-19;
const double mu0 = 0.003162278;
const double M_m = 3.0e-26; // Massa molecul
const double MH2O = 1.8e-2;
const double De = 2.3e-11; // 2.3e-11;
const double Pi = 3.1415;
const double sigma = 0.0; // 0.072;
const double Na = 6.0e23; // ����� ��������

const double Tm_0 = 1.0e6;

double P, P_i, P_f, Pg, dP, P_eq, Rho_g, P_corr;
double tau;
double Tm, T_i, T_out, T_out_r, T_d;
double T_010, T_005, T_001;
double Rb, RRb, RRRb, Mg;
double Nu, Nu_0;
double Cp_i, Cp_f, Cp, dCp, Cp_int, Cp_eq;
double E, eps, beta, D_eff;
double V_0, N_0, J, W_cr, V_m, d_m, J_e, Nb;

double *CpR, *NuR, *r, *CpR_n;
int i, Im, l, out_num, out_num_r;
double dr, r_tmp;

int Nu_flag, Cp_flag, Nucl_flag, Relaugh_flag;

char out_name[25];
FILE *in_file, *out_file, *Cp_file, *cut_file, *T_file;
char s[128];


// Non-dimensionalization

double t_nd, R_nd, dp_i, dC_i, gamma, epsilon, Nu_nd, m_nd;