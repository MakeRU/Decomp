// Decomp.cpp : 

#include <iostream>
#include <cstdio>
#include <math.h>
#include <time.h>


#include "Data.h"

double D_dihot(double eps)
{
	double D_min, D_max, D_int, D_tmp, dXi, Xi, Eq_tmp;
	int N, i;

	N = 10000000;
	D_min = 0.0;
	D_max = 1.0e10;
	dXi = 1.0e-4;

	D_tmp = (D_max + D_min) / 2.0;
	D_tmp = D_max;
	do {
		D_int = 0.0;
		for (i = 0; i < N; i++)
		{
			Xi = 1.0 + i * dXi;
			D_int = D_int + (exp(-(Xi*Xi / 2.0 + 1.0 / Xi - 3.0 / 2.0)*D_tmp / 2.0) / (Xi*Xi))*dXi;
		}
		Eq_tmp = D_tmp / 2.0*D_int - eps;
		if (Eq_tmp < 0.0) { D_min = D_tmp; }
		else { D_max = D_tmp; };
		D_tmp = (D_max + D_min) / 2.0;
	} while (D_max - D_min > 1.0e-10);

	return D_tmp;
}

int main()
{

	Nucl_flag = 1;

	in_file = fopen( "Init.txt", "rt" );

	fgets(s, 128, in_file);
	fgets(s, 128, in_file); sscanf(s, "%lf", &Tm);
	fgets(s, 128, in_file); sscanf(s, "%lf", &P_i);
	fgets(s, 128, in_file); sscanf(s, "%lf", &T_d);
	fgets(s, 128, in_file); sscanf(s, "%lf", &T_out);
	fgets(s, 128, in_file); sscanf(s, "%lf", &T_out_r);
	fgets(s, 128, in_file); sscanf(s, "%lf", &tau);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Im);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Cp_flag);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Nu_flag);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Relaugh_flag);
	
	fclose(in_file);

	T_i = 0.0;

	T_010 = Tm;
	T_005 = Tm;
	T_001 = Tm;


	P_i = P_i*p0;
	P_f = 1.0*p0;
	P = P_i;

	if (Nucl_flag == 1) {
		Rb = 1.0e-8;
		RRb = 0.0;
		Nb = 1.0e12;
		Pg = P_i + 2.0*sigma / Rb;
		Mg = Pg * 4.0*Pi *  Rb * Rb * Rb * MH2O / (3.0* Rg * T );
	}

	Cp_i = Kh * sqrt(P_i);
	Cp_f = Kh * sqrt(P_f);
	Cp = Kh * sqrt(Pg);

	P_corr = P_i;
	P = P_corr;

//	Nu_0 = 1.0e1;
	E = E0 * (1 - 12.0 * Cp_i);
	Nu_0 = mu0 * exp(E / (kB*T));
		
	out_num = 0;
	out_num_r = 0;

	// Non-dimensionalization

	dp_i = P_i - P_f;
	dC_i = Cp_i - Cp_f;
	Nu_nd = Nu_0;
	t_nd = Nu_nd / dp_i;
	R_nd = sqrt(t_nd*De);
	m_nd = 4.0 * Pi / 3.0 * R_nd * R_nd * R_nd* P_i* MH2O / (Rg * T);


	Im = 50;
	dr = 49.0*Rb / Im;
	r = new double[Im+1];
	CpR = new double[Im+1];
	CpR_n = new double[Im+1];
	NuR = new double[Im+1];
	for (i = 0; i <= Im; i++)
	{
		r[i] = Rb + i * dr;
		CpR[i] = Cp_i;
		CpR_n[i] = CpR[i];
		E = E0 * (1 - 12.0 * CpR[i]);
		NuR[i] = mu0 * exp(E / (kB*T));
	}


	out_file = fopen("Data/ViscBubble.dat", "wt");
	fprintf(out_file, "# P_i = %10.8lf \n", P_i / p0);
	fprintf(out_file, "# P_f = %10.8lf \n", P_f / p0);

	do {
//		printf("Time %3.8lf s \n", T_i);

		if (T_i > out_num * T_out) {
			fprintf(out_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf  \t %10.8lf \t %e \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \n",
				T_i, P / p0, (Pg - 2.0*sigma / Rb )/ p0, Rb*1.0e6, RRb*1.0e6, Mg, Nu, Cp, J, Nb);
			out_num = out_num + 1;
		}


		if (T_i < T_d) { P = P_corr - T_i * (P_corr - P_f) / (T_d); }
		else { P = P_f; };
		Cp_eq = Kh * sqrt(P);

		if (Nucl_flag == 0)
		{
			P_eq = Cp * Cp / (Kh*Kh);
			dP = P_eq - P;
			W_cr = 16.0*Pi*sigma*sigma*sigma / (3.0*dP*dP);
			V_0 = 1.0e-15;
			N_0 = rho0 * Cp / M_m;
			V_m = 0.018 / (1000 * Na);
			d_m = pow(6.0 / (Pi * N_0), 1.0 / 3.0);
			J_e = 2.0*N_0*N_0*V_m*De / d_m * sqrt(sigma / (kB*T));
			J = J_e * exp(-W_cr / (kB*T));
			Nb = Nb + V_0 * J*tau;
			if (Nb > 1.0) {
				Nb = Nb / V_0;
				Rb = 2 * sigma / (dP);
				Pg = P_i;
				Mg = 4.0*Pi * Pg * Rb * Rb * Rb * MH2O / (3.0* Rg * T);
				Rho_g = Pg * MH2O / (Rg * T);
				Cp = Cp_i - Nb * Mg / rho0;
				if (Cp < 0) { Cp = Kh * sqrt(P); };


				eps = rho0 * (Cp_i - Cp) / Rho_g;
				if (eps > 10.0) { D_eff = 12.0*De*eps*eps / Pi; }
				else { D_eff = 2.0*De*eps; };

				D_eff = De * D_dihot(eps);


				for (i = 0; i <= Im; i++) {
					r[i] = Rb + i * dr;
				}

				for (i = 0; i <= Im; i++) {
					CpR[i] = Cp;
					E = E0 * (1 - 12.0 * CpR[i]);
					NuR[i] = mu0 * exp(E / (kB*T));
				}

				if (Nu_flag == 0) { Nu = Nu_0; }
				else { Nu = NuR[0]; };

				RRb = (Pg - P)*Rb / (4.0*Nu);

				Nucl_flag = 1;
			}
		}

		if (Nucl_flag == 1)
		{
			dr = 49.0*Rb / Im;
			for (i = 0; i <= Im; i++)
			{
				r[i] = Rb + i * dr;
			}

			CpR[0] = Kh * sqrt(Pg);
			CpR[Im] = Cp_i;
			for (i = 1; i < Im; i++)
			{
			// CpR_n[i] = CpR[i]  + ( (2.0*De /r[i] - (RRb*Rb*Rb)/(r[i]*r[i]))*(CpR[i+1]-CpR[i-1])/(2.0*dr) +
			//			De *(CpR[i+1]-2.0*CpR[i]+ CpR[i-1])/(dr*dr) )*tau;
			//	CpR_n[i] = CpR[0] + (r[i]-r[0]) *(CpR[Im]-CpR[0])/(r[Im]-r[0]);
				CpR_n[i] = CpR[i]  + ( (2.0*De /r[i] - (RRb*Rb*Rb)/(r[i]*r[i]))*(CpR[i+1]-CpR[i])/(dr) +
										De *(CpR[i+1]-2.0*CpR[i]+ CpR[i-1])/(dr*dr) )*tau;
			}
			CpR_n[0] = CpR[0];
			CpR_n[Im] = Cp_i;

			for (i = 0; i <= Im; i++)
			{
			CpR[i] = CpR_n[i];
			}

			for (i = 0; i <= Im; i++)
			{
				E = E0 * (1 - 12.0 * CpR[i]);
				NuR[i] = mu0 * exp(E / (kB*T));
			}

			Nu = 0.0;
			for (i = 0; i < Im; i++)
			{
				r_tmp = (r[i + 1] + r[i]) / 2.0;
				Nu = Nu + (NuR[i+1] + NuR[i]) / 2.0 * dr / (r_tmp * r_tmp * r_tmp * r_tmp);
			}

			if (Nu_flag == 0) { Nu = Nu_0; }
			else { Nu = 3.0* Rb*Rb*Rb*Nu; }

			if (Cp_flag == 0) { dCp = Cp - Cp_eq; }
			else { dCp = (CpR[1] - CpR[0]) / dr; };
//			if (dCp < 0.0) { dCp = 0.0; }
//			if (Pg > P_eq)  { dCp = 0.0; }


			if (Relaugh_flag == 0)
			{
			Mg = Mg + 4.0*Pi* Rb * Rb * rho0 * De * dCp *tau;
			Pg = P;
			Rb=pow(3.0*Mg*Rg*T/(4.0*Pi*Pg*MH2O), 1.0/3.0);
			Rho_g = Pg * MH2O / (Rg * T);
			Cp = Cp_i - Nb * Mg / rho0;
			Cp_eq = Kh * sqrt(P);
			if (Cp < Cp_eq) { Cp = Kh * sqrt(P); };				
			} else 
			{
			RRb = (Pg - P - 2.0*sigma / Rb)*Rb / (4.0*Nu);
//			RRRb = (Pg - P - 2.0*sigma/Rb - 4.0*Nu*RRb/Rb) / (Rb*rho0) - 3.0*RRb*RRb/(2.0*Rb);
//			RRb = RRb + RRRb * tau;
//			if (RRb < 0.0) {RRb = 0.0;};
			Rb = Rb + RRb * tau;
			Mg = Mg + 4.0*Pi* Rb * Rb * rho0 * De * dCp *tau;
			Pg = Mg * 3.0* Rg * T / (4.0*Pi *  Rb * Rb * Rb * MH2O);
			Rho_g = Pg * MH2O / (Rg * T);
			Cp = Cp_i - Nb * Mg / rho0;
			Cp_eq = Kh * sqrt(P);
			if (Cp < Cp_eq) { Cp = Kh * sqrt(P); };
			}

			if (T_i > out_num_r * T_out_r) {
				if (out_num_r < 10000000) { sprintf(out_name, "Data/%d.dat", out_num_r); };
				if (out_num_r < 1000000) { sprintf(out_name, "Data/0%d.dat", out_num_r); };
				if (out_num_r < 100000) { sprintf(out_name, "Data/00%d.dat", out_num_r); };
				if (out_num_r < 10000) { sprintf(out_name, "Data/000%d.dat", out_num_r); };
				if (out_num_r < 1000) { sprintf(out_name, "Data/0000%d.dat", out_num_r); };
				if (out_num_r < 100) { sprintf(out_name, "Data/00000%d.dat", out_num_r); };
				if (out_num_r < 10) { sprintf(out_name, "Data/000000%d.dat", out_num_r); };
				cut_file = fopen(out_name, "wt");
				for (i = 0; i <= Im; i++) {
					fprintf(cut_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \n", r[i] * 1.0e6, r[i] / r[0], CpR[i], NuR[i]);
					//	fprintf(Cp_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \n", T_i, r[i] * 1.0e6, r[i]/r[0], CpR[i], NuR[i]);
				}
				fclose(cut_file);
				out_num_r = out_num_r + 1;
				printf("Time %3.8lf s \n", T_i);
			}

		}


		T_i = T_i + tau;

	} while (T_i <= Tm + tau);

	fclose(out_file);

	return 0;
}

