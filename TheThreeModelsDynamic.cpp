#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <array>
#include <complex>
#include <windows.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

using namespace std;

typedef complex<double> doublec;

#define npart 1
#define neq (2*npart)

//#define nstep 1000

constexpr int numRuns = 1;

#define LLG_DYN 1
//#undef LLG_DYN
#define LLG_TIMING 1
#undef LLG_TIMING
#define SW_DYN 1
#undef SW_DYN
#define SW_TIMING 1
#undef SW_TIMING
#define SW_APPROX 1
#undef SW_APPROX
#define SW_APPROX_TIMING 1
#undef SW_APPROX_TIMING

constexpr int n_max_vec = 30;

constexpr double Pi = 3.1415926535897932384626433832795;
constexpr double Pix2 = 6.283185307179586476925286766559;
constexpr double Pis2 = 1.5707963267948966192313216916398;

constexpr double uround = 1.0e-8;

struct sAnizo {
public: double ax, ay, az;
};
struct sReadData {
public: double x, y, z, theta_ea, phi_ea, volume, k, Msat, theta_sol, phi_sol;
};
struct sCoef {
public: int vecin;
	  double C, coef, xx, xy, xz, yx, yy, yz, zx, zy, zz;
};
struct Camp {
public: double H, theta, phi, Hx, Hy, Hz;
};
struct Moment {
public: double M, theta, phi, Mx, My, Mz;
};
struct SW_function_params {
public: double Hk_fcn, H_fcn, theta0_fcn;
};

struct Camp H_de_t(double t);
double SafeAcos(double x);
void position_coeficients(struct sReadData D1, struct sReadData D2, struct sCoef* Pos_Coef, double* dist);
void function_neighbours(void);
void anisotropy_coef(void);
int fcn(double t, const double y[], double yprime[], void* params);
//void fcn(int n_equ, double t, double y[], double yprime[]);
int fcn_xyz(double t, const double y[], double yprime[], void* params);
//void fcn_xyz(int n_equ, double t, double y[], double yprime[]);
double stability_test(double solutie[], double solutie_old[]);
double SW_fcn(double x, void* params);
void SW_angle_single_full(int part, double* sol, struct sReadData Medium_single, struct Camp H);
void SW_approx(int part, double* sol, struct sReadData Medium_single, struct Camp H);
double interpolare_GA(double A, double theta_in, double theta_tar);
void timp_2D_3D(double t, double h, double* sol_old, double* sol_target, double* sol_new);

//******************************************************

constexpr double VolumTotal = 1 * npart;

constexpr double alpha = 0.01;
constexpr double miu0 = 4.0e-7 * M_PI;

constexpr double Ms = 795774.7154594767; // 2.668e5;
constexpr double K1 = 1.0e5;

const double gamma = 2.210173e5;
const double time_norm = (1.0 + alpha * alpha) / gamma / Ms;

double theta_ea = 0.1 * M_PI / 180.0;
double phi_ea = 0.1 * M_PI / 180.0;

double theta_h = 20.0 * M_PI / 180.0;
double phi_h = 0.1 * M_PI / 180.0;

double theta_0 = theta_ea + 1.0e-4;
double phi_0 = phi_ea + 1.0e-4;

const double Hk = 2.0 * fabs(K1) / miu0 / Ms / Ms;

constexpr double	Field_period = 100000.0;					// 1e5 picosec -> 1e-7 s		1e4 picosec -> 1e-8 s
constexpr double	Freq_H = 1.0 / Field_period;			// 1e5 ps -> 1e7 Hz (10 MHz)	1e4 ps -> 1e8 Hz (100 MHz)
constexpr int		nperiods = 1;
constexpr double	t_max = nperiods * Field_period;

constexpr double	Read_period = Field_period / 1000.0;	// 1e1 picosec -> 1e-11 s
constexpr double	Freq_Read = 1.0 / Read_period;			// 1e1 ps -> 1e11 Hz (100 GHz)

//vrem sa citim nperiods of Field_period cu Read_period
constexpr int		nstep = nperiods * (int)(Field_period / Read_period);
constexpr double	t_step = t_max / (nstep - 1);		// reading time step

double T_Larmor = 1.0 / (gamma * 2.0 * fabs(K1) / (miu0 * Ms) / (2.0 * M_PI)) * 1.0e12;

const double Exch = 0.1;
const double prag_vecini = 5.0e-3;

//******************************************************

static std::array<struct sAnizo, npart>Anizo_params{};
static std::array<struct sReadData, npart>Medium{};
static std::array<struct Camp, nstep> H{};
static struct Camp Hext;
static std::array<std::array<struct Moment, npart>, nstep> M{};
struct Moment Msys;

static std::array<double, npart> Hx_part{};
static std::array<double, npart> Hy_part{};
static std::array<double, npart> Hz_part{};

static std::array<int, npart> neighbours{};
static std::array<std::array<struct sCoef, n_max_vec>, npart> Position_Coef{};

//******************************************************

char save_file_1_LLG[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2020\\SW---LLG\\Timing\\LLG_time_Js1-K1e5_th20_10MHz-MHL.dat";
char save_file_2_SW[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2020\\SW---LLG\\Timing\\SW_time_Js1-K1e5_th20_100MHz-MHL.dat";
char save_file_3_SWAPPROX[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2020\\SW---LLG\\Timing\\SWAPPROX_time_Js1-K1e5_th20_100MHz-MHL.dat";

//const double FieldMax = +1193662.0 / Ms;
//const double FieldMin = -1193662.0 / Ms;

const double FieldMax = +3.0 * Hk;
const double FieldMin = -3.0 * Hk;

//******************************************************

int main()
{
	int i, j;

	double       t = 0.0;
	//double       tend = 2000.0;
	//double		 tstep = 1.0 / (time_norm * 1.0e12);
	double		 last_t = 0.0;
	double y[neq];
	double y_old[neq];
	double y_target[neq];

	for (i = 0; i < nstep; i++)
	{
		H[i].H = FieldMax * cos(i * nperiods * 2.0 * Pi / (nstep - 1));
		H[i].theta = theta_h;
		H[i].phi = phi_h;
		H[i].Hx = H[i].H * sin(H[i].theta) * cos(H[i].phi);
		H[i].Hy = H[i].H * sin(H[i].theta) * sin(H[i].phi);
		H[i].Hz = H[i].H * cos(H[i].theta);
	}

	for (i = 0; i < npart; i++)
	{
		Medium[i].x = 0.0;
		Medium[i].y = 0.0;
		Medium[i].z = 0.0;
		Medium[i].volume = 1.0;
		Medium[i].Msat = Ms;
		Medium[i].k = 1.0;
		Medium[i].theta_ea = theta_ea;
		Medium[i].phi_ea = phi_ea;
		Medium[i].theta_sol = theta_0;
		Medium[i].phi_sol = phi_0;
	}

	for (i = 0; i < npart; i++)
	{
		y[2 * i + 0] = Medium[i].theta_sol;
		y[2 * i + 1] = Medium[i].phi_sol;
	}

	anisotropy_coef();
	function_neighbours();

	//////////////////////////////////////////////////////////////////////////
	// MHL LLG
	//////////////////////////////////////////////////////////////////////////
#ifdef LLG_DYN
	{
		FILE* fp;
		int num_LLG = 0;

		fp = fopen(save_file_1_LLG, "w");
		fclose(fp);

		gsl_odeiv2_system sys = { fcn, NULL, neq, NULL };

		double ug_th = theta_h;
		double ug_ph = phi_h;
		double MHL_projection;

		//initial conditions
		for (i = 0; i < npart; i++)
		{
			y[2 * i + 0] = theta_h;
			y[2 * i + 1] = phi_h;
		}

		t = 0.0;
		last_t = 0.0;

		for (int h = 0; h < nstep; h++)
		{
			for (j = 0; j < neq; j++)
				y_old[j] = y[j];

			Hext = H[h];
			num_LLG = 0;

			gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-8, 1e-8);

			fp = fopen(save_file_1_LLG, "a");


			{
				gsl_odeiv2_driver_apply(d, &t, t + t_step, y);

				Msys.Mx = 0.0; Msys.My = 0.0; Msys.Mz = 0.0;

				for (j = 0; j < npart; j++)
				{
					Msys.Mx += Medium[j].volume * sin(y[2 * j + 0]) * cos(y[2 * j + 1]);
					Msys.My += Medium[j].volume * sin(y[2 * j + 0]) * sin(y[2 * j + 1]);
					Msys.Mz += Medium[j].volume * cos(y[2 * j + 0]);
				}

				num_LLG++;

				MHL_projection = (Msys.Mx * sin(Hext.theta) * cos(Hext.phi) + Msys.My * sin(Hext.theta) * sin(Hext.phi) + Msys.Mz * cos(Hext.theta));
			}

			fprintf(fp, "%20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf\n", t * time_norm * 1.0e12, Hext.Hx, Hext.Hy, Hext.Hz, Msys.Mx, Msys.My, Msys.Mz, Hext.H, MHL_projection);

			last_t = (t * time_norm * 1.0e12);

			fclose(fp);
			gsl_odeiv2_driver_free(d);

			printf("%07.4lf ps -> H = %07.4lf -> M = %07.4lf \n", t * time_norm * 1.0e12, Hext.H, MHL_projection);
		}
	}
#endif
#ifdef LLG_TIMING
	{
		DWORD starttime, elapsedtime;
		starttime = timeGetTime();
		for (int numMHLs = 0; numMHLs < numRuns; numMHLs++)
		{
			int num_LLG = 0;

			gsl_odeiv2_system sys = { fcn, NULL, neq, NULL };

			double ug_th = theta_h;
			double ug_ph = phi_h;
			double MHL_projection;

			//initial conditions
			for (i = 0; i < npart; i++)
			{
				y[2 * i + 0] = theta_h;
				y[2 * i + 1] = phi_h;
			}

			t = 0.0;
			last_t = 0.0;

			for (int h = 0; h < nstep; h++)
			{
				for (j = 0; j < neq; j++)
					y_old[j] = y[j];

				Hext = H[h];
				num_LLG = 0;

				gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-8, 1e-8);

				while (((stability_test(y, y_old)) > 1.0e-4) || (num_LLG < 100))
				{
					gsl_odeiv2_driver_apply(d, &t, t + tstep, y);

					Msys.Mx = 0.0; Msys.My = 0.0; Msys.Mz = 0.0;

					for (j = 0; j < npart; j++)
					{
						Msys.Mx += Medium[j].volume * sin(y[2 * j + 0]) * cos(y[2 * j + 1]);
						Msys.My += Medium[j].volume * sin(y[2 * j + 0]) * sin(y[2 * j + 1]);
						Msys.Mz += Medium[j].volume * cos(y[2 * j + 0]);
					}

					num_LLG++;

					MHL_projection = (Msys.Mx * sin(Hext.theta) * cos(Hext.phi) + Msys.My * sin(Hext.theta) * sin(Hext.phi) + Msys.Mz * cos(Hext.theta));
				}

				last_t = (t * time_norm * 1.0e12);

				gsl_odeiv2_driver_free(d);
			}
		}
		elapsedtime = timeGetTime() - starttime;
		printf("Time Elapsed %10d mSecs \n", (int)elapsedtime);
	}
#endif
	//////////////////////////////////////////////////////////////////////////
	// MHL SW
	//////////////////////////////////////////////////////////////////////////
#ifdef SW_DYN
	{
		FILE* fp;

		fp = fopen(save_file_2_SW, "w");
		fclose(fp);

		double ug_th = theta_h;
		double ug_ph = phi_h;
		double MHL_projection;

		//initial conditions
		//////////////////////////////////// SATURATE!
		for (i = 0; i < npart; i++)
		{
			Medium[i].theta_sol = H[0].theta;
			Medium[i].phi_sol = H[0].phi;
			y[2 * i + 0] = Medium[i].theta_sol;
			y[2 * i + 1] = Medium[i].phi_sol;
			y_old[2 * i + 0] = y[2 * i + 0];
			y_old[2 * i + 1] = y[2 * i + 1];
		}

		Hext = H[0];
		for (j = 0; j < npart; j++)
			SW_angle_single_full(j, &y[2 * j + 0], Medium[j], Hext);

		for (j = 0; j < npart; j++)
		{
			Medium[j].theta_sol = y[2 * j + 0];
			Medium[j].phi_sol = y[2 * j + 1];
			y_old[2 * j + 0] = y[2 * j + 0];
			y_old[2 * j + 1] = y[2 * j + 1];
		}
		/////////////////////////////////////////////////

		for (int h = 0; h < nstep; h++)
		{
			Hext = H[h];
			if (H[h].H < 0.0)
			{
				printf("a");
			}

			for (j = 0; j < npart; j++)
			{
				y[2 * j + 0] = Medium[j].theta_sol;
				y[2 * j + 1] = Medium[j].phi_sol;
				y_old[2 * j + 0] = y[2 * j + 0];
				y_old[2 * j + 1] = y[2 * j + 1];
			}

			fp = fopen(save_file_2_SW, "a");

			for (j = 0; j < npart; j++)
				SW_angle_single_full(j, &y[2 * j + 0], Medium[j], Hext);

			/////////////////////////////////////////////////////////////
			for (j = 0; j < npart; j++)
			{
				y_target[2 * j + 0] = y[2 * j + 0];
				y_target[2 * j + 1] = y[2 * j + 1];
				y[2 * j + 0] = y_old[2 * j + 0];
				y[2 * j + 1] = y_old[2 * j + 1];
			}

			timp_2D_3D(t_step, Hext.H, y_old, y_target, y);
			///////////////////////////////////////////////////////////////

			Msys.Mx = 0.0; Msys.My = 0.0; Msys.Mz = 0.0;
			for (j = 0; j < npart; j++)
			{
				Medium[j].theta_sol = y[2 * j + 0];
				Medium[j].phi_sol = y[2 * j + 1];
				Msys.Mx += Medium[j].volume * sin(y[2 * j + 0]) * cos(y[2 * j + 1]);
				Msys.My += Medium[j].volume * sin(y[2 * j + 0]) * sin(y[2 * j + 1]);
				Msys.Mz += Medium[j].volume * cos(y[2 * j + 0]);
			}
			Msys.Mx /= VolumTotal;
			Msys.My /= VolumTotal;
			Msys.Mz /= VolumTotal;

			MHL_projection = (Msys.Mx * sin(Hext.theta) * cos(Hext.phi) + Msys.My * sin(Hext.theta) * sin(Hext.phi) + Msys.Mz * cos(Hext.theta));

			fprintf(fp, "%20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf\n", t * time_norm * 1.0e12, Hext.Hx, Hext.Hy, Hext.Hz, Msys.Mx, Msys.My, Msys.Mz, Hext.H, MHL_projection);

			fclose(fp);

			printf("%07.4lf ps -> H = %07.4lf -> M = %07.4lf \n", t * time_norm * 1.0e12, Hext.H, MHL_projection);
		}
	}
#endif
#ifdef SW_TIMING
	{
		DWORD starttime, elapsedtime;
		starttime = timeGetTime();
		for (int numMHLs = 0; numMHLs < numRuns; numMHLs++)
		{
			double ug_th = theta_h;
			double ug_ph = phi_h;
			double MHL_projection;

			//initial conditions
		//////////////////////////////////// SATURARTE!
			for (i = 0; i < npart; i++)
			{
				Medium[i].theta_sol = H[0].theta;
				Medium[i].phi_sol = H[0].phi;
				y[2 * i + 0] = Medium[i].theta_sol;
				y[2 * i + 1] = Medium[i].phi_sol;
				y_old[2 * i + 0] = y[2 * i + 0];
				y_old[2 * i + 1] = y[2 * i + 1];
			}

			Hext = H[0];
			for (j = 0; j < npart; j++)
				SW_angle_single_full(j, &y[2 * j + 0], Medium[j], Hext);

			for (j = 0; j < npart; j++)
			{
				Medium[j].theta_sol = y[2 * j + 0];
				Medium[j].phi_sol = y[2 * j + 1];
				y_old[2 * j + 0] = y[2 * j + 0];
				y_old[2 * j + 1] = y[2 * j + 1];
			}
			/////////////////////////////////////////////////

			for (int h = 0; h < nstep; h++)
			{
				Hext = H[h];

				for (j = 0; j < npart; j++)
				{
					y[2 * j + 0] = Medium[j].theta_sol;
					y[2 * j + 1] = Medium[j].phi_sol;
					y_old[2 * j + 0] = y[2 * j + 0];
					y_old[2 * j + 1] = y[2 * j + 1];
				}

				for (j = 0; j < npart; j++)
					SW_angle_single_full(j, &y[2 * j + 0], Medium[j], Hext);

				Msys.Mx = 0.0; Msys.My = 0.0; Msys.Mz = 0.0;
				for (j = 0; j < npart; j++)
				{
					Medium[j].theta_sol = y[2 * j + 0];
					Medium[j].phi_sol = y[2 * j + 1];
					Msys.Mx += Medium[j].volume * sin(y[2 * j + 0]) * cos(y[2 * j + 1]);
					Msys.My += Medium[j].volume * sin(y[2 * j + 0]) * sin(y[2 * j + 1]);
					Msys.Mz += Medium[j].volume * cos(y[2 * j + 0]);
				}
				Msys.Mx /= VolumTotal;
				Msys.My /= VolumTotal;
				Msys.Mz /= VolumTotal;

				MHL_projection = (Msys.Mx * sin(Hext.theta) * cos(Hext.phi) + Msys.My * sin(Hext.theta) * sin(Hext.phi) + Msys.Mz * cos(Hext.theta));
			}
		}
		elapsedtime = timeGetTime() - starttime;
		printf("Time Elapsed %10d mSecs \n", (int)elapsedtime);
	}
#endif
	//////////////////////////////////////////////////////////////////////////
	// MHL SW APPROX
	//////////////////////////////////////////////////////////////////////////
#ifdef SW_APPROX
	{
		FILE* fp;

		fp = fopen(save_file_3_SWAPPROX, "w");
		fclose(fp);

		double ug_th = theta_h;
		double ug_ph = phi_h;
		double MHL_projection;

		//initial conditions
	//////////////////////////////////// SATURARTE!
		for (i = 0; i < npart; i++)
		{
			Medium[i].theta_sol = H[0].theta;
			Medium[i].phi_sol = H[0].phi;
			y[2 * i + 0] = Medium[i].theta_sol;
			y[2 * i + 1] = Medium[i].phi_sol;
			y_old[2 * i + 0] = y[2 * i + 0];
			y_old[2 * i + 1] = y[2 * i + 1];
		}

		Hext = H[0];
		for (j = 0; j < npart; j++)
			SW_approx(j, &y[2 * j + 0], Medium[j], Hext);

		for (j = 0; j < npart; j++)
		{
			Medium[j].theta_sol = y[2 * j + 0];
			Medium[j].phi_sol = y[2 * j + 1];
			y_old[2 * j + 0] = y[2 * j + 0];
			y_old[2 * j + 1] = y[2 * j + 1];
		}
		/////////////////////////////////////////////////

		for (int h = 0; h < nstep; h++)
		{
			Hext = H[h];

			for (j = 0; j < npart; j++)
			{
				y[2 * j + 0] = Medium[j].theta_sol;
				y[2 * j + 1] = Medium[j].phi_sol;
				y_old[2 * j + 0] = y[2 * j + 0];
				y_old[2 * j + 1] = y[2 * j + 1];
			}

			fp = fopen(save_file_3_SWAPPROX, "a");

			for (j = 0; j < npart; j++)
				SW_approx(j, &y[2 * j + 0], Medium[j], Hext);

			Msys.Mx = 0.0; Msys.My = 0.0; Msys.Mz = 0.0;
			for (j = 0; j < npart; j++)
			{
				Medium[j].theta_sol = y[2 * j + 0];
				Medium[j].phi_sol = y[2 * j + 1];
				Msys.Mx += Medium[j].volume * sin(y[2 * j + 0]) * cos(y[2 * j + 1]);
				Msys.My += Medium[j].volume * sin(y[2 * j + 0]) * sin(y[2 * j + 1]);
				Msys.Mz += Medium[j].volume * cos(y[2 * j + 0]);
			}
			Msys.Mx /= VolumTotal;
			Msys.My /= VolumTotal;
			Msys.Mz /= VolumTotal;

			MHL_projection = (Msys.Mx * sin(Hext.theta) * cos(Hext.phi) + Msys.My * sin(Hext.theta) * sin(Hext.phi) + Msys.Mz * cos(Hext.theta));

			fprintf(fp, "%20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf\n", t * time_norm * 1.0e12, Hext.Hx, Hext.Hy, Hext.Hz, Msys.Mx, Msys.My, Msys.Mz, Hext.H, MHL_projection);

			fclose(fp);

			printf("%07.4lf ps -> H = %07.4lf -> M = %07.4lf \n", t * time_norm * 1.0e12, Hext.H, MHL_projection);
		}
	}
#endif
#ifdef SW_APPROX_TIMING
	{
		DWORD starttime, elapsedtime;
		starttime = timeGetTime();
		for (int numMHLs = 0; numMHLs < numRuns; numMHLs++)
		{
			double ug_th = theta_h;
			double ug_ph = phi_h;
			double MHL_projection;

			//initial conditions
		//////////////////////////////////// SATURARTE!
			for (i = 0; i < npart; i++)
			{
				Medium[i].theta_sol = H[0].theta;
				Medium[i].phi_sol = H[0].phi;
				y[2 * i + 0] = Medium[i].theta_sol;
				y[2 * i + 1] = Medium[i].phi_sol;
				y_old[2 * i + 0] = y[2 * i + 0];
				y_old[2 * i + 1] = y[2 * i + 1];
			}

			Hext = H[0];
			for (j = 0; j < npart; j++)
				SW_approx(j, &y[2 * j + 0], Medium[j], Hext);

			for (j = 0; j < npart; j++)
			{
				Medium[j].theta_sol = y[2 * j + 0];
				Medium[j].phi_sol = y[2 * j + 1];
				y_old[2 * j + 0] = y[2 * j + 0];
				y_old[2 * j + 1] = y[2 * j + 1];
			}
			/////////////////////////////////////////////////

			for (int h = 0; h < nstep; h++)
			{
				Hext = H[h];

				for (j = 0; j < npart; j++)
				{
					y[2 * j + 0] = Medium[j].theta_sol;
					y[2 * j + 1] = Medium[j].phi_sol;
					y_old[2 * j + 0] = y[2 * j + 0];
					y_old[2 * j + 1] = y[2 * j + 1];
				}

				for (j = 0; j < npart; j++)
					SW_approx(j, &y[2 * j + 0], Medium[j], Hext);

				Msys.Mx = 0.0; Msys.My = 0.0; Msys.Mz = 0.0;
				for (j = 0; j < npart; j++)
				{
					Medium[j].theta_sol = y[2 * j + 0];
					Medium[j].phi_sol = y[2 * j + 1];
					Msys.Mx += Medium[j].volume * sin(y[2 * j + 0]) * cos(y[2 * j + 1]);
					Msys.My += Medium[j].volume * sin(y[2 * j + 0]) * sin(y[2 * j + 1]);
					Msys.Mz += Medium[j].volume * cos(y[2 * j + 0]);
				}
				Msys.Mx /= VolumTotal;
				Msys.My /= VolumTotal;
				Msys.Mz /= VolumTotal;

				MHL_projection = (Msys.Mx * sin(Hext.theta) * cos(Hext.phi) + Msys.My * sin(Hext.theta) * sin(Hext.phi) + Msys.Mz * cos(Hext.theta));
			}
		}
		elapsedtime = timeGetTime() - starttime;
		printf("Time Elapsed %10d mSecs \n", (int)elapsedtime);
	}
#endif
	return(0);
}

//**************************************************************************
//**************************************************************************
//**************************************************************************
//**************************************************************************
//**************************************************************************

struct Camp H_de_t(double t)
{
	struct Camp H;

	H.theta = theta_h;
	H.phi = phi_h;

	double ampl = FieldMax;

	H.H = ampl * cos(Pix2 * t / Field_period);

	H.Hx = H.H * sin(H.theta) * cos(H.phi);
	H.Hy = H.H * sin(H.theta) * sin(H.phi);
	H.Hz = H.H * cos(H.theta);

	return H;
}

//**************************************************************************

double SafeAcos(double x)
{
	if (x < -1.0) x = -1.0;
	else if (x > 1.0) x = 1.0;
	return acos(x);
}

//**************************************************************************

void position_coeficients(struct sReadData D1, struct sReadData D2, sCoef* Pos_Coef, double* dist)
{
	*dist = sqrt((D2.x - D1.x) * (D2.x - D1.x) + (D2.y - D1.y) * (D2.y - D1.y) + (D2.z - D1.z) * (D2.z - D1.z));
	double r = (double)rand() / RAND_MAX;

	if (*dist == 0)
	{
		Pos_Coef->C = r * 0.1;

		Pos_Coef->coef = 0.0;
		Pos_Coef->xx = 0.0;
		Pos_Coef->xy = 0.0;
		Pos_Coef->xz = 0.0;
		Pos_Coef->yx = 0.0;
		Pos_Coef->yy = 0.0;
		Pos_Coef->yz = 0.0;
		Pos_Coef->zx = 0.0;
		Pos_Coef->zy = 0.0;
		Pos_Coef->zz = 0.0;
	}
	else
	{
		double rx = (D2.x - D1.x) / *dist;
		double ry = (D2.y - D1.y) / *dist;
		double rz = (D2.z - D1.z) / *dist;

		Pos_Coef->C = 0.0;

		Pos_Coef->coef = D2.volume * D2.Msat / 4.0 / M_PI / *dist / *dist / *dist / Ms;

		Pos_Coef->xx = 3.0 * rx * rx - 1.0;
		Pos_Coef->xy = 3.0 * rx * ry;
		Pos_Coef->xz = 3.0 * rx * rz;
		Pos_Coef->yx = Pos_Coef->xy;
		Pos_Coef->yy = 3.0 * ry * ry - 1.0;
		Pos_Coef->yz = 3.0 * ry * rz;
		Pos_Coef->zx = Pos_Coef->xz;
		Pos_Coef->zy = Pos_Coef->yz;
		Pos_Coef->zz = 3.0 * rz * rz - 1.0;
	}
}

//**************************************************************************

void function_neighbours(void)
{
	sReadData Data1, Data2;
	sCoef Pos_Coef;
	double distance;
	int neighbours_max = 0, neighbours_med = 0;
	int f;

	for (int i = 0; i < npart; i++)
	{
		neighbours[i] = 0;
		Data1 = Medium[i];

		for (f = 0; f < npart; f++)
		{
			Data2 = Medium[f];

			if ((i != f))
			{
				position_coeficients(Data1, Data2, &Pos_Coef, &distance);
				if (distance < 1.0e-9)
				{
					if (neighbours[i] + 1 > n_max_vec - 1)
					{
						printf("\n\n PREA MULTI VECINI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n\n");
						getchar();
						break;
					}
					else
					{
						neighbours[i]++;
						Pos_Coef.vecin = f;  //<<<---- asta e vecin
						Position_Coef[i][neighbours[i] - 1] = Pos_Coef;
					}
					continue; // din acelasi cub
				}

				if (Pos_Coef.coef > prag_vecini)
				{
					if (neighbours[i] + 1 > n_max_vec - 1)
					{
						printf("\n\n PREA MULTI VECINI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n\n");
						getchar();
						break;
					}
					else
					{
						neighbours[i]++;
						Pos_Coef.vecin = f;  //<<<---- asta e vecin
						Position_Coef[i][neighbours[i] - 1] = Pos_Coef;
					}
				}
			}
		}

		printf("particula: %d  cu  %d  vecini \n", i, neighbours[i]);
		neighbours_med += neighbours[i];
		if (neighbours[i] > neighbours_max) neighbours_max = neighbours[i];
	}
	printf("Numar maxim de vecini: %d      Numar mediu de vecini: %f \n", neighbours_max, (float)neighbours_med / npart);
}

//**************************************************************************

void anisotropy_coef(void)
{
	for (int i = 0; i < npart; i++)
	{
		Anizo_params[i].ax = sin(Medium[i].theta_ea) * cos(Medium[i].phi_ea);
		Anizo_params[i].ay = sin(Medium[i].theta_ea) * sin(Medium[i].phi_ea);
		Anizo_params[i].az = cos(Medium[i].theta_ea);
	}
}

//**************************************************************************

int fcn(double time, const double input[], double deriv[], void* params)
//void fcn(int n_equ, double t, double input[], double deriv[])
{
	int i, j, k, kk;
	double s0, c0, s1, c1;
	double term1, term2, term3, factor;
	double mp_x, mp_y, mp_z;

	for (i = 0; i < npart; i++)                 //  Fcn_Dopri(x,y,k1), sol->y era param in Dopri
	{
		k = 2 * i;

		s0 = sin(input[k]);
		c0 = cos(input[k]);
		s1 = sin(input[k + 1]);
		c1 = cos(input[k + 1]);

		term1 = Anizo_params[i].ax * s0 * c1 + Anizo_params[i].ay * s0 * s1 + Anizo_params[i].az * c0;
		term2 = -Anizo_params[i].ax * s0 * s1 + Anizo_params[i].ay * s0 * c1;
		term3 = Anizo_params[i].ax * c0 * c1 + Anizo_params[i].ay * c0 * s1 - Anizo_params[i].az * s0;
		factor = 1.0 /*+ 2.0*K2*term1*term1 / K1 + 3.0*K3*term1*term1*term1*term1 / K1*/;

		Hx_part[i] = Hext.Hx + Medium[i].k * Hk * factor * (-term1 * term2 * s1 / s0 + term1 * term3 * c0 * c1);
		Hy_part[i] = Hext.Hy + Medium[i].k * Hk * factor * (term1 * term2 * c1 / s0 + term1 * term3 * c0 * s1);
		Hz_part[i] = Hext.Hz + Medium[i].k * Hk * factor * (-term1 * term3 * s0);

		for (j = 0; j < neighbours[i]; j++)
		{
			kk = 2 * Position_Coef[i][j].vecin;
			mp_x = sin(input[kk]) * cos(input[kk + 1]);
			mp_y = sin(input[kk]) * sin(input[kk + 1]);
			mp_z = cos(input[kk]);

			Hx_part[i] += Position_Coef[i][j].C * Hk * mp_x + Position_Coef[i][j].coef * (Position_Coef[i][j].xx * mp_x + Position_Coef[i][j].xy * mp_y + Position_Coef[i][j].xz * mp_z);
			Hx_part[i] += Position_Coef[i][j].C * Hk * mp_y + Position_Coef[i][j].coef * (Position_Coef[i][j].yx * mp_x + Position_Coef[i][j].yy * mp_y + Position_Coef[i][j].yz * mp_z);
			Hx_part[i] += Position_Coef[i][j].C * Hk * mp_z + Position_Coef[i][j].coef * (Position_Coef[i][j].zx * mp_x + Position_Coef[i][j].zy * mp_y + Position_Coef[i][j].zz * mp_z);
		}

		deriv[k] = Hx_part[i] * (alpha * c0 * c1 - s1) + Hy_part[i] * (alpha * c0 * s1 + c1) - Hz_part[i] * alpha * s0;
		deriv[k + 1] = (-Hx_part[i] * (alpha * s1 + c0 * c1) + Hy_part[i] * (alpha * c1 - c0 * s1) + Hz_part[i] * s0) / s0;
	}

	return GSL_SUCCESS;
}

//**************************************************************************

int fcn_xyz(double time, const double input[], double deriv[], void* params)
//void fcn_xyz(int n_equ, double t, double input[], double deriv[])
{
	int i, j, k, kk;
	double aux;
	double mp_x, mp_y, mp_z;

	for (i = 0; i < npart; i++)                 //  Fcn_Dopri(x,y,k1), sol->y era param in Dopri
	{
		k = 3 * i;

		aux = Anizo_params[i].ax * input[k + 0] + Anizo_params[i].ay * input[k + 1] + Anizo_params[i].az * input[k + 2];

		Hx_part[i] = Hext.Hx - Medium[i].k * Hk * aux * Anizo_params[i].ax;
		Hy_part[i] = Hext.Hy - Medium[i].k * Hk * aux * Anizo_params[i].ax;
		Hz_part[i] = Hext.Hz - Medium[i].k * Hk * aux * Anizo_params[i].ax;

		for (j = 0; j < neighbours[i]; j++)
		{
			kk = 3 * Position_Coef[i][j].vecin;
			mp_x = input[kk + 0];
			mp_y = input[kk + 1];
			mp_z = input[kk + 2];

			Hx_part[i] += Position_Coef[i][j].C * Hk * mp_x + Position_Coef[i][j].coef * (Position_Coef[i][j].xx * mp_x + Position_Coef[i][j].xy * mp_y + Position_Coef[i][j].xz * mp_z);
			Hy_part[i] += Position_Coef[i][j].C * Hk * mp_y + Position_Coef[i][j].coef * (Position_Coef[i][j].yx * mp_x + Position_Coef[i][j].yy * mp_y + Position_Coef[i][j].yz * mp_z);
			Hz_part[i] += Position_Coef[i][j].C * Hk * mp_z + Position_Coef[i][j].coef * (Position_Coef[i][j].zx * mp_x + Position_Coef[i][j].zy * mp_y + Position_Coef[i][j].zz * mp_z);
		}
		//(Hy1*mz(t)                - Hz1*my(t))              - alpha1*(Hy1*mx(t)*my(t)                         + Hz1*mx(t)*mz(t)                        - Hx1*my(t)*my(t)                          - Hx1*mz(t)*mz(t))
		deriv[k + 0] = (Hy_part[i] * input[k + 2] - Hz_part[i] * input[k + 1]) - alpha * (Hy_part[i] * input[k + 0] * input[k + 1] + Hz_part[i] * input[k + 0] * input[k + 2] - Hx_part[i] * input[k + 1] * input[k + 1] - Hx_part[i] * input[k + 2] * input[k + 2]);

		//(Hz1*mx(t)             - Hx1*mz(t))                 - alpha1*(Hx1*mx(t)*my(t)                         + Hz1*my(t)*mz(t)                        - Hy1*mx(t)*mx(t)                          - Hy1*mz(t)*mz(t))
		deriv[k + 1] = (Hz_part[i] * input[k + 0] - Hx_part[i] * input[k + 2]) - alpha * (Hx_part[i] * input[k + 0] * input[k + 1] + Hz_part[i] * input[k + 1] * input[k + 2] - Hy_part[i] * input[k + 0] * input[k + 0] - Hy_part[i] * input[k + 2] * input[k + 2]);

		//(Hx1*my(t)              - Hy1*mx(t))                - alpha1*(Hx1*mx(t)*mz(t)                         + Hy1*my(t)*mz(t)                       - Hz1*mx(t)*mx(t)                           - Hz1*my(t)*my(t))
		deriv[k + 2] = (Hx_part[i] * input[k + 1] - Hy_part[i] * input[k + 0]) - alpha * (Hx_part[i] * input[k + 0] * input[k + 2] + Hy_part[i] * input[k + 1] * input[k + 2] - Hz_part[i] * input[k + 0] * input[k + 0] - Hz_part[i] * input[k + 1] * input[k + 1]);
	}

	return GSL_SUCCESS;
}

//**************************************************************************

double stability_test(double solutie[], double solutie_old[])
{
	double diference = 0.0;
	double torq, torq_old;
	double proj, proj_old;
	double th, ph, th_old, ph_old;

	for (int i = 0; i < npart; i++)
	{
		if (fabs(solutie_old[2 * i + 0] - solutie[2 * i + 0]) > diference)
			diference = fabs(solutie_old[2 * i + 0] - solutie[2 * i + 0]);

		if (fabs(solutie_old[2 * i + 1] - solutie[2 * i + 1]) > diference)
			diference = fabs(solutie_old[2 * i + 1] - solutie[2 * i + 1]);

		//solutie_old[2 * i + 0] = solutie[2 * i + 0];
		//solutie_old[2 * i + 1] = solutie[2 * i + 1];

		th = atan(tan(solutie[2 * i + 0]));
		ph = atan(tan(solutie[2 * i + 1]));

		torq = (th * Hext.phi - ph * Hext.theta) * (th * Hext.phi - ph * Hext.theta) + (ph - Hext.phi) * (ph - Hext.phi) + (Hext.theta - th) * (Hext.theta - th);
		proj = (sin(th) * cos(ph) * sin(Hext.theta) * cos(Hext.phi) + sin(th) * sin(ph) * sin(Hext.theta) * sin(Hext.phi) + cos(th) * cos(Hext.theta));

		th_old = atan(tan(solutie_old[2 * i + 0]));
		ph_old = atan(tan(solutie_old[2 * i + 1]));

		//solutie_old[2 * i + 0] = solutie[2 * i + 0];
		//solutie_old[2 * i + 1] = solutie[2 * i + 1];

		torq_old = (th_old * Hext.phi - ph_old * Hext.theta) * (th_old * Hext.phi - ph_old * Hext.theta) + (ph_old - Hext.phi) * (ph_old - Hext.phi) + (Hext.theta - th_old) * (Hext.theta - th_old);
		proj_old = (sin(th_old) * cos(ph_old) * sin(Hext.theta) * cos(Hext.phi) + sin(th_old) * sin(ph_old) * sin(Hext.theta) * sin(Hext.phi) + cos(th_old) * cos(Hext.theta));

		if (fabs(torq_old - torq) > diference)
			diference = fabs(torq_old - torq);

		solutie_old[2 * i + 0] = solutie[2 * i + 0];
		solutie_old[2 * i + 1] = solutie[2 * i + 1];

		//  		if (fabs(proj_old - proj) > diference)
		//  			diference = fabs(proj_old - proj);
	}

	return diference;
}

//**************************************************************************

double SW_fcn(double x, void* params)
{
	struct SW_function_params* p = (struct SW_function_params*) params;

	double Hk_fcn = p->Hk_fcn;
	double H_fcn = p->H_fcn;
	double theta0_fcn = p->theta0_fcn;

	return Hk_fcn * cos(x) * sin(x) + H_fcn * sin(x - theta0_fcn);
}

//**************************************************************************

void SW_angle_single_full(int part, double* sol, struct sReadData Medium_single, struct Camp H)
{
	static double loco_angle_M, loco_old_angle_M;
	static double loco_angle_H;
	static double Hc;
	int n_found = 0, state = 1;
	double g, temp_angle_H, temp_Hea, temp_Hha;
	double hx, hy, hz;
	double mx, my, mz;

	int status, n_sols, ram = 0;
	int iter = 0, max_iter = 100;

	//double* sols;
	//double* limits;
	std::array<double, 4> sols{};
	std::array<double, 8> limits{};

	const gsl_root_fsolver_type* T;
	gsl_root_fsolver* s;
	double r = 0;
	double x_lo, x_hi;
	gsl_function F;

	H.H = sqrt(H.Hx * H.Hx + H.Hy * H.Hy + H.Hz * H.Hz);
	H.theta = SafeAcos(H.Hz / H.H);
	H.phi = atan2(H.Hy, H.Hx);

	hx = sin(H.theta) * cos(H.phi);
	hy = sin(H.theta) * sin(H.phi);
	hz = cos(H.theta);

	//	Convert_sol_to_loco(sol, loco_angle_M);
	{
		mx = sin(sol[0]) * cos(sol[1]);
		my = sin(sol[0]) * sin(sol[1]);
		mz = cos(sol[0]);
		loco_angle_M = SafeAcos(Anizo_params[part].ax * mx + Anizo_params[part].ay * my + Anizo_params[part].az * mz);
	}

	loco_angle_H = SafeAcos(Anizo_params[part].ax * hx + Anizo_params[part].ay * hy + Anizo_params[part].az * hz);
	double its_a_sin = sin(loco_angle_H);
	double its_a_cos = cos(loco_angle_H);
	double aux1 = pow(its_a_sin * its_a_sin, (1.0 / 3.0));
	double aux2 = pow(its_a_cos * its_a_cos, (1.0 / 3.0));
	g = pow(aux1 + aux2, -1.5);
	Hc = Medium_single.k * Hk * g;

	//parametrii functiei de rezolvat (in coordonate loco)
	double Hk_fcn = Medium_single.k * Hk;
	double H_fcn = H.H;
	double theta0_fcn = loco_angle_H;
	struct SW_function_params params = { Medium_single.k * Hk, H.H, loco_angle_H };
	////////////////////////////////////////////////////
	F.function = &SW_fcn;
	F.params = &params;
	T = gsl_root_fsolver_brent;
	/////////////////////////////////////////////////////////

	loco_old_angle_M = SafeAcos(Anizo_params[part].ax * sin(sol[0]) * cos(sol[1]) + Anizo_params[part].ay * sin(sol[0]) * sin(sol[1]) + Anizo_params[part].az * cos(sol[0]));
	if (loco_old_angle_M < Pis2)
	{
		state = 1;
	}
	else
	{
		state = -1;
	}

	if (H.H == 0.0)
	{
		//sols = (double*)calloc(1, sizeof(double));
		//limits = (double*)calloc(1, sizeof(double));
		if (state > 0)
		{
			loco_angle_M = 0.0;
		}
		else
		{
			loco_angle_M = Pi;
		}
	}
	else
	{
		temp_Hea = H.H * cos(loco_angle_H);
		temp_Hha = H.H * sin(loco_angle_H);
		if ((temp_Hea >= 0.0) && (temp_Hha >= 0.0))
		{
			temp_angle_H = Pi + atan(-pow(tan(loco_angle_H) * tan(loco_angle_H), (1.0 / 6.0)));	//C1
		}
		if ((temp_Hea >= 0.0) && (temp_Hha < 0.0))
		{
			temp_angle_H = Pi + atan(pow(tan(loco_angle_H) * tan(loco_angle_H), (1.0 / 6.0)));	//C4
		}
		if ((temp_Hea < 0.0) && (temp_Hha < 0.0))
		{
			temp_angle_H = Pix2 + atan(-pow(tan(loco_angle_H) * tan(loco_angle_H), (1.0 / 6.0)));	//C3
		}
		if ((temp_Hea < 0.0) && (temp_Hha >= 0.0))
		{
			temp_angle_H = atan(pow(tan(loco_angle_H) * tan(loco_angle_H), (1.0 / 6.0)));	//C2
		}

		if (H.H >= Hc)
			n_sols = 2;
		else
			n_sols = 4;

		//sols = (double*)calloc(n_sols, sizeof(double));
		//limits = (double*)calloc(2 * n_sols, sizeof(double));

		//printf("%d: temp_angle: %lf      (Hea, Hha):  (%lf, %lf)\n", n_sols, temp_angle_H, temp_Hea, temp_Hha);

		if (fabs(H.H) >= Hc)
		{
			if ((temp_Hea >= 0.0) && (temp_Hha >= 0.0))
			{
				limits[0] = 0.0; limits[1] = Pis2; limits[2] = Pi; limits[3] = 3.0 * Pis2;	//C1
			}
			if ((temp_Hea >= 0.0) && (temp_Hha < 0.0))
			{
				limits[0] = Pis2; limits[1] = Pi; limits[2] = 3.0 * Pis2; limits[3] = Pix2;	//C4
			}
			if ((temp_Hea < 0.0) && (temp_Hha < 0.0))
			{
				limits[0] = 0.0; limits[1] = Pis2; limits[2] = Pi; limits[3] = 3.0 * Pis2;	//C3
			}
			if ((temp_Hea < 0.0) && (temp_Hha >= 0.0))
			{
				limits[0] = Pis2; limits[1] = Pi; limits[2] = 3.0 * Pis2; limits[3] = Pix2;	//C2
			}
		}
		else
		{
			if ((temp_Hea >= 0.0) && (temp_Hha >= 0.0))
			{
				limits[0] = 0.0; limits[1] = Pis2; limits[2] = Pi; limits[3] = 3.0 * Pis2; limits[4] = Pis2; limits[5] = temp_angle_H; limits[6] = temp_angle_H; limits[7] = Pi;
			}
			if ((temp_Hea >= 0.0) && (temp_Hha < 0.0))
			{
				limits[0] = Pis2; limits[1] = Pi; limits[2] = 3.0 * Pis2; limits[3] = Pix2; limits[4] = Pi; limits[5] = temp_angle_H; limits[6] = temp_angle_H; limits[7] = 3.0 * Pis2;
			}
			if ((temp_Hea < 0.0) && (temp_Hha < 0.0))
			{
				limits[0] = 0.0; limits[1] = Pis2; limits[2] = Pi; limits[3] = 3.0 * Pis2; limits[4] = 3.0 * Pis2; limits[5] = temp_angle_H; limits[6] = temp_angle_H; limits[7] = Pix2;
			}
			if ((temp_Hea < 0.0) && (temp_Hha >= 0.0))
			{
				limits[0] = Pis2; limits[1] = Pi; limits[2] = 3.0 * Pis2; limits[3] = Pix2; limits[4] = 0.0; limits[5] = temp_angle_H; limits[6] = temp_angle_H; limits[7] = Pis2;
			}
		}

		for (int scontor = 0; scontor < n_sols; scontor++)
		{
			s = gsl_root_fsolver_alloc(T);
			x_lo = limits[2 * scontor];
			x_hi = limits[2 * scontor + 1];
			gsl_root_fsolver_set(s, &F, x_lo, x_hi);
			do
			{
				iter++;
				status = gsl_root_fsolver_iterate(s);
				r = gsl_root_fsolver_root(s);
				x_lo = gsl_root_fsolver_x_lower(s);
				x_hi = gsl_root_fsolver_x_upper(s);
				status = gsl_root_test_interval(x_lo, x_hi, 0, 1.0e-7);
			} while (status == GSL_CONTINUE && iter < max_iter);
			sols[scontor] = r;
			gsl_root_fsolver_free(s);
		}

		if (H.H >= Hc)
		{
			for (int scontor = 0; scontor < n_sols; scontor++)
			{
				if ((Hk_fcn * cos(2.0 * sols[scontor]) + H_fcn * cos(sols[scontor] - theta0_fcn)) > 0)
					loco_angle_M = sols[scontor];
			}
		}
		else
		{
			int contor_sol_bune = 0;
			//double sol_bune[2];
			std::array<double, 2> sol_bune{};
			for (int scontor = 0; scontor < n_sols; scontor++)
			{
				if ((Hk_fcn * cos(2.0 * sols[scontor]) + H_fcn * cos(sols[scontor] - theta0_fcn)) > 0)
				{
					sol_bune[contor_sol_bune] = sols[scontor];
					contor_sol_bune++;
				}
			}

			if (contor_sol_bune != 2)
			{
				printf("SOMETHING FISHY!\n");
			}

			if (fabs(loco_old_angle_M - sol_bune[0]) > fabs(loco_old_angle_M - sol_bune[1]))
				loco_angle_M = sol_bune[1];
			else
				loco_angle_M = sol_bune[0];
		}
	}

	//	Convert_loco_to_sol(loco_angle_M, sol, H);
	{
		double /*hx, hy, hz, mx, my, mz, */ug_ea_H, temp_ea, temp_H;
		hx = sin(H.theta) * cos(H.phi);
		hy = sin(H.theta) * sin(H.phi);
		hz = cos(H.theta);
		if (H.H < 0.0)
		{
			hx *= (-1.0);
			hy *= (-1.0);
			hz *= (-1.0);
		}

		ug_ea_H = SafeAcos(Anizo_params[part].ax * hx + Anizo_params[part].ay * hy + Anizo_params[part].az * hz);

		if (loco_angle_M <= Pi)
		{
			temp_ea = sin(ug_ea_H - loco_angle_M) / sin(ug_ea_H);
			temp_H = sin(loco_angle_M) / sin(ug_ea_H);
		}
		else
		{
			temp_ea = sin(ug_ea_H - (Pix2 - loco_angle_M)) / sin(ug_ea_H);
			temp_H = sin(Pix2 - loco_angle_M) / sin(ug_ea_H);
		}

		mx = (temp_ea * Anizo_params[part].ax + temp_H * hx);
		my = (temp_ea * Anizo_params[part].ay + temp_H * hy);
		mz = (temp_ea * Anizo_params[part].az + temp_H * hz);

		sol[0] = SafeAcos(mz);

		if ((mx == 0.0) && (my == 0.0))
		{
			sol[1] = 0.0;
		}
		else
		{
			sol[1] = SafeAcos(mx / sqrt(mx * mx + my * my));
			if (my < 0)
			{
				sol[1] = Pix2 - sol[1];
			}
		}
	}
	// 	free(sols);
	// 	free(limits);
}

//**************************************************************************

void SW_approx(int part, double* sol, struct sReadData Medium_single, struct Camp H)//double th_h, double ph_h, double a, )
{
	double g;
	double its_a_sin, its_a_cos, aux1, aux2;
	double hx, hy, hz, mx, my, mz;
	double ug_ea_H, ug_ea_M, temp_ea, temp_H;

	double loco_angle_H;
	double loco_old_angle_M;
	double loco_angle_M;
	double altr_angle_M;
	double Hc;
	double sol_bune[2];

	//camp(th_h, ph_h, a, &H);
	//Add_interactions(part, &H, sol);
	//H.Hamp = (H.Hamp < 1.0e-4) ? 1.1e-4 : H.Hamp;
	//amplitudini[part] = H.Hamp;

	H.H = sqrt(H.Hx * H.Hx + H.Hy * H.Hy + H.Hz * H.Hz);
	H.theta = SafeAcos(H.Hz / H.H);
	H.phi = atan2(H.Hy, H.Hx);

	hx = sin(H.theta) * cos(H.phi);
	hy = sin(H.theta) * sin(H.phi);
	hz = cos(H.theta);
	loco_angle_H = SafeAcos(Anizo_params[part].ax * hx + Anizo_params[part].ay * hy + Anizo_params[part].az * hz);


	//	Convert_sol_to_loco(sol, loco_angle_M);
	mx = sin(sol[0]) * cos(sol[1]);
	my = sin(sol[0]) * sin(sol[1]);
	mz = cos(sol[0]);
	loco_old_angle_M = SafeAcos(Anizo_params[part].ax * mx + Anizo_params[part].ay * my + Anizo_params[part].az * mz);


	its_a_sin = sin(loco_angle_H);
	its_a_cos = cos(loco_angle_H);
	aux1 = pow(its_a_sin * its_a_sin, (1.0 / 3.0));
	aux2 = pow(its_a_cos * its_a_cos, (1.0 / 3.0));
	g = pow(aux1 + aux2, -1.5);
	Hc = Medium_single.k * Hk * g;

	double Hk_fcn = Medium_single.k * Hk;
	double H_fcn = H.H;
	double h = H_fcn / Hk_fcn;
	double phi = loco_angle_H;

	hy = h * cos(phi);
	hx = h * sin(phi);

	double d = 1 - h * h;

	doublec e(0.0, 0.0), fp(0.0, 0.0), fm(0.0, 0.0), mp(0.0, 0.0), mm(0.0, 0.0), mpx(0.0, 0.0), mmx(0.0, 0.0), t1(0.0, 0.0), t2(0.0, 0.0);

	e = d * cos(acos(doublec(54.0 * hx * hx * hy * hy / d / d / d - 1.0)) / 3.0);
	fp = sqrt(9.0 * hy * hy + 6.0 * d + 6.0 * e);
	fm = -fp;
	mp = (fp + sqrt(2.0 * fp * fp - 18.0 * e + 54.0 * hy * (1.0 + hx * hx) / fp)) / 6.0 - hy / 2.0;
	mm = (fm - sqrt(2.0 * fm * fm - 18.0 * e + 54.0 * hy * (1.0 + hx * hx) / fm)) / 6.0 - hy / 2.0;
	mpx = sqrt(1.0 - mp * mp);
	mmx = sqrt(1.0 - mm * mm);
	t1 = atan(mpx / mp);
	t2 = atan(mmx / mm);
	sol_bune[0] = real(t1);
	sol_bune[1] = Pi + real(t2);

	if (H.H == 0.0) {
		loco_angle_M = (loco_old_angle_M < Pis2) ? 0.0 : Pi;
	}
	if (H.H >= Hc) {
		loco_angle_M = (phi < Pis2) ? sol_bune[0] : sol_bune[1];
	}
	else {
		loco_angle_M = (fabs(loco_old_angle_M - sol_bune[0]) > fabs(loco_old_angle_M - sol_bune[1])) ? sol_bune[1] : sol_bune[0];
		altr_angle_M = (fabs(loco_old_angle_M - sol_bune[0]) > fabs(loco_old_angle_M - sol_bune[1])) ? sol_bune[0] : sol_bune[1];
	}

	// 	if (dezint == 1) {
	// 		for (int part = 0; part < npart; part++) {
	// 			if (H_SW[part].Hamp < Hc[part]) {
	// 				Eng_barM = energy_barrier_wood(part, loco_angle_M, altr_angle_M, loco_angle_H, H_SW, Hc);
	// 				tau[part] = 1.0e-12 * P[part].T_Larmor_Hk * exp(Eng_barM / kBT);
	// 				prob = 1.0 - exp(-t_univ / tau[part]);
	// 				r = uniform_rand(0.0, 1.0);
	// 				if (r < prob)
	// 					loco_angle_M[part] = altr_angle_M[part];
	// 			}
	// 		}
	// 	}

	ug_ea_H = loco_angle_H;
	ug_ea_M = loco_angle_M;
	temp_ea = (loco_angle_M <= Pi) ? sin(ug_ea_H - ug_ea_M) / sin(ug_ea_H) : sin(ug_ea_H - (Pix2 - ug_ea_M)) / sin(ug_ea_H);
	temp_H = (loco_angle_M <= Pi) ? sin(ug_ea_M) / sin(ug_ea_H) : sin(Pix2 - ug_ea_M) / sin(ug_ea_H);

	hx = sin(H.theta) * cos(H.phi);
	hy = sin(H.theta) * sin(H.phi);
	hz = cos(H.theta);
	mx = (temp_ea * Anizo_params[part].ax + temp_H * hx);
	my = (temp_ea * Anizo_params[part].ay + temp_H * hy);
	mz = (temp_ea * Anizo_params[part].az + temp_H * hz);

	sol[0] = SafeAcos(mz);
	sol[1] = atan2(my, mx);
}

//**************************************************************************

double interpolare_GA(double A, double theta_in, double theta_tar)
{
	/*if (fabs(theta_in - theta_tar) < 1.0e-5)
		return 20.0;*/

	double p1, p2, p3, p4, Y0, ampl, expo;
	double th_in = 0, theta_in_c, tmc, u, x, timp_polynom = 0;

	p1 = 1.58519 - 0.08438 * pow(0.64213, A);
	p2 = (1.18986 + 1.85715 * pow(0.62603, A)) * theta_tar;
	p3 = (-0.17513 - 1.72152 * pow(0.62553, A)) * theta_tar * theta_tar;
	p4 = (0.03716 + 0.36513 * pow(0.62556, A)) * theta_tar * theta_tar * theta_tar;

	Y0 = 16.91794 + 90.13333 * pow(0.71009, A);
	ampl = 3.85951 + 196.57648 * pow(0.48572, A);
	expo = exp(-0.5 * (theta_tar - Pis2) * (theta_tar - Pis2) / 0.64 / 0.64);

	theta_in_c = p1 + p2 + p3 + p4;
	u = 2 * theta_in_c - theta_tar;
	tmc = Y0 + ampl * expo;

	if (theta_in < theta_tar)
		th_in = 2 * theta_tar - theta_in;
	if (theta_in > u)
		th_in = 2 * u - theta_in;
	if (theta_in > theta_tar&& theta_in < u)
		th_in = theta_in;

	x = th_in - theta_in_c;
	timp_polynom = (x > 0) ? 0.7417 + 0.906 * exp(x / 0.5193) : -(0.7417 + 0.906 * exp(-x / 0.5193));

	return (timp_polynom + tmc);
}

//**************************************************************************

void timp_2D_3D(double t, double h, double *sol_old, double *sol_target, double *sol_new)
{
	double dphi;
	double raport;
	double mx0, mx1, my0, my1, mz0, mz1, mxf, myf, mzf;
	double mx0r, my0r, mz0r, mx1r, my1r, mz1r;
	double t_max, loco_angle_with_ea_old, loco_angle_with_ea_target, angle_target, angle_new;
	double th_m0r, th_m1r, ph_m0r, ph_m1r;
	double s0, c0, s1, c1;
	double loco_angle_old_uthea, loco_angle_target_uthea;

	t_max = 0;

	for (int i = 0; i < npart; i++) {

		mx0 = sin(sol_old[2 * i + 0]) * cos(sol_old[2 * i + 1]);
		my0 = sin(sol_old[2 * i + 0]) * sin(sol_old[2 * i + 1]);
		mz0 = cos(sol_old[2 * i + 0]);

		mx1 = sin(sol_target[2 * i + 0]) * cos(sol_target[2 * i + 1]);
		my1 = sin(sol_target[2 * i + 0]) * sin(sol_target[2 * i + 1]);
		mz1 = cos(sol_target[2 * i + 0]);

		loco_angle_with_ea_old = SafeAcos(Anizo_params[i].ax * mx0 + Anizo_params[i].ay * my0 + Anizo_params[i].az * mz0);
		angle_target = SafeAcos(mx0 * mx1 + my0 * my1 + mz0 * mz1);
		loco_angle_with_ea_target = SafeAcos(Anizo_params[i].ax * mx1 + Anizo_params[i].ay * my1 + Anizo_params[i].az * mz1);
		double utheax = cos(Medium[i].theta_ea) * cos(Medium[i].phi_ea);
		double utheay = cos(Medium[i].theta_ea) * sin(Medium[i].phi_ea);
		double utheaz = -sin(Medium[i].theta_ea);
		loco_angle_old_uthea = SafeAcos(mx0 * utheax + my0 * utheay + mz0 * utheaz);
		loco_angle_target_uthea = SafeAcos(mx1 * utheax + my1 * utheay + mz1 * utheaz);
		if ((loco_angle_target_uthea<Pis2 && loco_angle_old_uthea>Pis2) || (loco_angle_target_uthea > Pis2&& loco_angle_old_uthea < Pis2))
			loco_angle_with_ea_old = Pix2 - loco_angle_with_ea_old;


		// se roteste pozitia veche a lui H astfel incat acesta sa fie dupa Oz si se afla noul M
		double R1[3][3], R2[3][3];
		
		s0 = sin(sol_target[2 * i + 0]);
		s1 = sin(sol_target[2 * i + 1]);
		c0 = cos(sol_target[2 * i + 0]);
		c1 = cos(sol_target[2 * i + 1]);

		R1[0][0] = c0 * c1;  R1[0][1] = c0 * s1;  R1[0][2] = -s0;
		R1[1][0] = -s1;       R1[1][1] = c1;       R1[1][2] = 0.0;
		R1[2][0] = s0 * c1;  R1[2][1] = s0 * s1;  R1[2][2] = c0;

		R2[0][0] = R1[0][0]; R2[0][1] = R1[1][0]; R2[0][2] = R1[2][0];
		R2[1][0] = R1[0][1]; R2[1][1] = R1[1][1]; R2[1][2] = R1[2][1];
		R2[2][0] = R1[0][2]; R2[2][1] = R1[1][2]; R2[2][2] = R1[2][2];

		mx0r = R1[0][0] * mx0 + R1[0][1] * my0 + R1[0][2] * mz0;
		my0r = R1[1][0] * mx0 + R1[1][1] * my0 + R1[1][2] * mz0;
		mz0r = R1[2][0] * mx0 + R1[2][1] * my0 + R1[2][2] * mz0;

		th_m0r = SafeAcos(mz0r);
		ph_m0r = atan2(my0r, mx0r);

		// se modifica th_m0r si ph_m0r conform regulilor

		double T_Larmor_Hk = Pix2 / (gamma * Hk * Medium[i].k * Ms) * 1.0e+12;

		raport = FieldMax/*fabs(h)*/ / Hk;
		raport = (raport < 2.0) ? 2.0 : raport;

		dphi = 0.46 + 0.976 * raport;
		ph_m1r = ph_m0r + dphi * t * Pix2 / T_Larmor_Hk;

		t_max = T_Larmor_Hk * interpolare_GA(raport, loco_angle_with_ea_old, loco_angle_with_ea_target) / 4.0;

		angle_new = (t < t_max) ? angle_target * 4.0 * (1.0 + fabs(cos(angle_target))) * t / t_max : angle_target;
		th_m1r = th_m0r - angle_new;
		if (th_m1r < 0.0)
			th_m1r = 0.0;

		// se roteste totul inapoi rezultand noua solutie
		mx1r = sin(th_m1r) * cos(ph_m1r);
		my1r = sin(th_m1r) * sin(ph_m1r);
		mz1r = cos(th_m1r);

		mxf = R2[0][0] * mx1r + R2[0][1] * my1r + R2[0][2] * mz1r;
		myf = R2[1][0] * mx1r + R2[1][1] * my1r + R2[1][2] * mz1r;
		mzf = R2[2][0] * mx1r + R2[2][1] * my1r + R2[2][2] * mz1r;

		sol_new[2 * i + 0] = SafeAcos(mzf);
		sol_new[2 * i + 1] = atan2(myf, mxf);
	}

}

//**************************************************************************

