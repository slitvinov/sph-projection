#ifndef ZS_PARTICLES_H
#define ZS_PARTICLES_H
#include "ZS_Physics.h"

class ZS_Particles : public ZS_Physics
{
	public:

		vector< vector<double> >	particles_r;
		vector< vector<double> >	particles_v;
		vector< vector<double> >	particles_a;		
		vector<double>				particles_rho;
		vector<double>				particles_p;
		vector<double>				particles_m;

		vector<double>				particles_w;
		vector< vector<double> >	particles_r_projection;
		vector< vector<double> >	particles_v_projection;
		vector< vector<double> >	particles_a_projection;		
		vector<double>				particles_rho_projection;
		vector<double>				particles_p_projection;
		vector<double>				particles_m_projection;

		long						particles_num;
		
	public:

		ZS_Particles();
		~ZS_Particles();

		void particles_init();
		void particles_projection_position();
		void particles_read(FILE *fid);
		void particles_write(FILE *fid);

};

#endif //ZS_PARTICLES_H
