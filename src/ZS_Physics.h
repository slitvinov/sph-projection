#ifndef ZS_PHYSICS_H
#define ZS_PHYSICS_H
#include "ZS_Control.h"
#include <vector>

class ZS_Physics : public ZS_Control
{
	public:

		int					    phy_num_dim;               // number of dimension
		long 					phy_step_start;            // start of step
		long 					phy_step_end;              // end of step
		long                    phy_msd_step;
		double                  phy_overlap;               // overlap for cutoff radience
		double					phy_cut_off;               // cutoff radience of normal particles
		double					phy_cut_off_macro;         // cutoff radience of macro-particles
		double					phy_dt;                    // time step of user's setting
		double                  phy_dt_CFL;                // time step corresponding with CFL condition
		double                  phy_tau_dt;                // factor of Kolmogorov time
		double                  phy_tau_macro;             // relaxation time in stochastic model
		double					phy_t;                     // simulation time
		double					phy_time_start;            // start of simulation time
		double					phy_time_end;              // end of simulation time
		double					phy_v0;					   // initial velocity for normal-particles (or micro-particles)
		double					phy_v0_macro;			   // initial velocity for macro-particles
		double					phy_rho;                   
		double					phy_eta;                   // viscosity
		double					phy_ksai;                  // bulk viscosity
		double					phy_temp;                  
		double					phy_k0;                    // initial kinetic temperature of micro-particles
		double					phy_c;                     // speed of sound
		double					phy_rho_0;                 // initial density for normal-particles (or micro-particles)
		double					phy_rho_0_macro;		   // initial density for macro-particles
		double					phy_b;                     // parameter in equation of state
		double					phy_gamma;				   // parameter in equation of state
		double					phy_gamma_dot;             // shear rate for lees-edwards boundary
		double					phy_Fc;                    // compressible force in forced turbulence
		double					phy_Fr;                    // rotational force in forced turbulence
		double                  phy_epsilon;               // energy dissipation rate
		double                  phy_epsilon_0;             // initial(average) energy dissipation rate
		double                  phy_A;                     // factor of deterministic force
		double                  phy_u_rms;                 // mean square root velocity of micro-particles
		double                  phy_u_rms_macro;           // mean square root velocity of macro-particles
		double                  phy_kolm_tau;              // Kolmogorov's time
                double                  phy_kolm_eta;              // Kolmogorov's scale
		double                  phy_Ma;                    // Mach number
		long                    phy_random_seed;
		double                  phy_kinetic_energy;
		  double                  phy_zeta_max;
		  double                 phy_zeta_min;
		vector<double>			phy_num_particles;         // number of normal particles in each direction
		vector<double>			phy_num_macroparticles;    // number of macro-particles in each direction
		vector<double>			phy_min_length;            // minimal value of box's length
		vector<double>			phy_max_length;            // maximal value of box's length
		vector<double>			phy_domain_length;         // the length of box
		vector<double>		        phy_body_force;            
		vector<double>                  phy_offset;
		vector<double>			phy_bcdef;
		vector<double>          phy_momentum;

	public:

		ZS_Physics();
		~ZS_Physics();

		void phy_init();
		void phy_read();
		void phy_print(FILE *fid);
	    void phy_get_vector(string cvector,vector<double> *vector_value);

		void get_time_step();
		void get_cut_off();
		void get_domain_length();
};

#endif //ZS_PHYSICS_H
