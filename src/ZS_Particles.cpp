#include "ZS_Particles.h"

ZS_Particles::ZS_Particles()
{
}

ZS_Particles::~ZS_Particles()
{
}

void ZS_Particles::particles_init()
{
	long i(0);
	double m(1);
	
	particles_num = 1;

	for(i=0;i < phy_num_dim;i++) 
	{
		particles_num =  particles_num * phy_num_particles[i];
		m = m * phy_domain_length[i] / phy_num_particles[i];
	}

	particles_r.resize(particles_num);
	particles_r_projection.resize(particles_num);
	
	particles_v.resize(particles_num);
	particles_v_projection.resize(particles_num);

	particles_a.resize(particles_num);
	particles_a_projection.resize(particles_num);
	
	for(i=0;i<particles_num;i++)
	{
		particles_r[i].resize(phy_num_dim,0.0);
		particles_v[i].resize(phy_num_dim,0.0);
		particles_a[i].resize(phy_num_dim,0.0);

		particles_v_projection[i].resize(phy_num_dim,0.0);
		particles_a_projection[i].resize(phy_num_dim,0.0);
	}

	particles_rho.resize(particles_num,0);
	particles_rho_projection.resize(particles_num,0);

	particles_m.resize(particles_num,m);
	particles_m_projection.resize(particles_num,m);

	particles_w.resize(particles_num, 0.0);
}

void ZS_Particles::particles_projection_position()
{
	int i(0),j(0),k(0),ndim(0);
	long id(0);
	vector<double> dx(phy_num_dim,0.0); 

	for(ndim=0;ndim<phy_num_dim;ndim++)
		dx[ndim] = phy_domain_length[ndim] / (phy_num_particles[ndim]);
		//dx[ndim] = phy_domain_length[ndim] / phy_num_particles[ndim];
	
	if (phy_num_dim == 2)
	{
		for(i=0;i<phy_num_particles[0];i++)
		{
			for(j=0;j<phy_num_particles[1];j++)
			{ 
				id = i * phy_num_particles[1] + j;
				if (ctr_projection_type == 0)
				{
					particles_r_projection[id].push_back(dx[0] / 2 + j * dx[0]);
					particles_r_projection[id].push_back(dx[1] / 2 + i * dx[1]);
					  
				}

				if (ctr_projection_type == 1)
				{
					particles_r_projection[id].push_back(particles_r[id][0]);
					particles_r_projection[id].push_back(particles_r[id][1]);
				}
			}
		}
	}

	if (phy_num_dim == 3)
	{
		for(i=0;i<phy_num_particles[0];i++)
		{
			for(j=0;j<phy_num_particles[1];j++)
			{ 
				for(k=0;k<phy_num_particles[2];k++)
				{
					id = i * phy_num_particles[1] * phy_num_particles[2] + j * phy_num_particles[1] + k;
					if (ctr_projection_type == 0)
					{
						particles_r_projection[id].push_back(dx[0] / 2 + k * dx[0]);
					        particles_r_projection[id].push_back(dx[1] / 2 + j * dx[1]);
					        particles_r_projection[id].push_back(dx[2] / 2 + i * dx[2]);
					      //cout << "particle " << id << " : Rx = " <<  particles_r_projection[id][0] << endl;
					      //cout << "particle " << id << " : Ry = " <<  particles_r_projection[id][1] << endl;
					      //cout << "particle " << id << " : Rz = " <<  particles_r_projection[id][2] << endl;
					     
					}

					if (ctr_projection_type == 1)
					{
						particles_r_projection[id].push_back(k * dx[0] + phy_offset[0] * dx[0]);
						particles_r_projection[id].push_back(j * dx[1] + phy_offset[1] * dx[1]);
					        particles_r_projection[id].push_back(i * dx[2] + phy_offset[2] * dx[2]);
					}
				}
			}
		}
	}
}

void ZS_Particles::particles_read(FILE *fid)
{
	long id(0);
	if (ctr_only_velocity == 1)
	{
		if (phy_num_dim == 2) 
		{
		
			double *valuePara = new double[4];	

			for(id=0;id < particles_num;id++)
			{
				fscanf(fid,"%lf %lf %lf %lf\n",&valuePara[0],&valuePara[1],&valuePara[2],&valuePara[3]);
				
				particles_r[id][0] = valuePara[0];
				particles_r[id][1] = valuePara[1];
				
				particles_v[id][0] = valuePara[2];
				particles_v[id][1] = valuePara[3];		
			}

			delete [] valuePara;
		}

		if (phy_num_dim == 3) 
		{
			double *valuePara = new double[6];

			for(id=0;id < particles_num;id++)
			{
				fscanf(fid,"%lf %lf %lf %lf %lf %lf\n",&valuePara[0],&valuePara[1],&valuePara[2],&valuePara[3],&valuePara[4],&valuePara[5]);
				
				particles_r[id][0] = valuePara[0];
				particles_r[id][1] = valuePara[1];
				particles_r[id][2] = valuePara[2];
				
				particles_v[id][0] = valuePara[3];
				particles_v[id][1] = valuePara[4];	
				particles_v[id][2] = valuePara[5];
				
			}

			delete [] valuePara;
		}
	}

	else
	{
		if (phy_num_dim == 2) 
		{
		
			double *valuePara = new double[7];	

			for(id=0;id < particles_num;id++)
			{
				fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf",&valuePara[0],&valuePara[1],&valuePara[2],&valuePara[3],&valuePara[4],&valuePara[5],&valuePara[6]);
				
				particles_r[id][0] = valuePara[0];
				particles_r[id][1] = valuePara[1];
				
				particles_v[id][0] = valuePara[2];
				particles_v[id][1] = valuePara[3];		
				
				particles_a[id][0] = valuePara[4];
				particles_a[id][1] = valuePara[5];

				particles_rho[id] = valuePara[6];
			}

			delete [] valuePara;
		}

		if (phy_num_dim == 3) 
		{
			double *valuePara = new double[10];

			for(id=0;id < particles_num;id++)
			{
				fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&valuePara[0],&valuePara[1],&valuePara[2],&valuePara[3],&valuePara[4],&valuePara[5],&valuePara[6],&valuePara[7],&valuePara[8],&valuePara[9]);
				
				particles_r[id][0] = valuePara[0];
				particles_r[id][1] = valuePara[1];
				particles_r[id][2] = valuePara[2];
				
				particles_v[id][0] = valuePara[3];
				particles_v[id][1] = valuePara[4];	
				particles_v[id][2] = valuePara[5];
				
				particles_a[id][0] = valuePara[6];
				particles_a[id][1] = valuePara[7];
				particles_a[id][2] = valuePara[8];

				particles_rho[id] = valuePara[9];
			}

			delete [] valuePara;
		}
	}
}

void ZS_Particles::particles_write(FILE *fid)
{
	long i(0);

	if (phy_num_dim == 2)
	{
		for(i=0;i<particles_num;i++)
			fprintf(fid,"%e %e %e %e %e %e %e\n",particles_r_projection[i][0],particles_r_projection[i][1],particles_v_projection[i][0],particles_v_projection[i][1],particles_a_projection[i][0],particles_a_projection[i][1],particles_rho_projection[i]);
		
		//fprintf(fid,"\n");
	}
	
	if (phy_num_dim == 3)
	{
		for(i=0;i<particles_num;i++)
			fprintf(fid,"%e %e %e %e %e %e %e %e %e %e\n",particles_r_projection[i][0],particles_r_projection[i][1],particles_r_projection[i][2],particles_v_projection[i][0],particles_v_projection[i][1],particles_v_projection[i][2],particles_a_projection[i][0],particles_a_projection[i][1],particles_a_projection[i][2],particles_rho_projection[i]);
	}

}
