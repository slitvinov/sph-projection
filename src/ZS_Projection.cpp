#include "ZS_Projection.h"

ZS_Projection::ZS_Projection()
{
}

ZS_Projection::~ZS_Projection()
{
}

void ZS_Projection::projection_init()
{
	ctrl_read();

        phy_init();
	
	particles_init();
	
	
}

void ZS_Projection::projection_calculate()
{
	long i(0);
	FILE *fid_in,*fid_out;
	  //errno_t err;
	char buffer[10];
	string input_name;
    
	i = phy_step_start;

	while( (i + 1)  <= phy_step_end)
	{
	  // cout << "n = " << i << ",    dt = " << phy_dt << endl;
		phy_t = i * phy_dt * ctr_output_freq;
		cout << "n = " << i << ",    t = " << phy_t << endl;

		input_name = ctr_input_name;
		sprintf(buffer,"%ld",i);
		input_name.append(buffer);

		  if( (fid_in  = fopen(input_name.c_str(), "rt" )) == NULL ) 
			cout << "Can't open input file !" << endl;

		ctr_output_name = input_name;
	
		ctr_output_name.append(".prj");

		
		if( (fid_out  = fopen(ctr_output_name.c_str(), "wt" )) == NULL ) 
			cout << "Can't open output file !" << endl;
		

		  //  cout << "ctr_only_velocity : " << ctr_only_velocity << endl;
		particles_read(fid_in);
		  // cout << "read!!" << endl;
		particles_projection_position();
		  //cout << "position!!" << endl;
		cell_list_init();
		  //cout << " cell init!!" << endl;
		boundary_init();
		cell_list_generate();
		boundary_generate();		
		projection_parameter();
		particles_write(fid_out);
		projection_parameter_reset();
		  //cout << "after reset: " << particles_rho_projection[16] <<endl;
		i++;
		fclose(fid_in);
		fclose(fid_out);
	}

	
}

void ZS_Projection::projection_parameter()
{
	if (ctr_test_case == 1)
		projection_parameter_pb();

	if (ctr_test_case == 2)
		projection_parameter_mb();

	if (ctr_test_case == 3)
		projection_parameter_pb();

	if (ctr_test_case == 4)
		projection_parameter_pb();

	if (ctr_test_case == 5)
		projection_parameter_pb();

	if (ctr_test_case == 6)
		projection_parameter_pb();
	
}

void ZS_Projection::projection_parameter_pb()
{
	if (ctr_list_type == 2)
	{
		if (phy_num_dim == 2)
		{
			for(cell_index[1]=0;cell_index[1]<cell_num[1];(cell_index[1])++)
				for(cell_index[0]=0;cell_index[0]<cell_num[0];(cell_index[0])++)
				{
					cell_id = cell_index[0] + cell_index[1] * cell_num[0];

					for(neighbour_index[1]=cell_index[1]-1;neighbour_index[1]<=cell_index[1]+1;(neighbour_index[1])++)
						for(neighbour_index[0]=cell_index[0]-1;neighbour_index[0]<=cell_index[0]+1;(neighbour_index[0])++)
						{
							boundary_EDS();
							
							neighbour_id = ((neighbour_index[0] + cell_num[0])%cell_num[0])  + ((neighbour_index[1] + cell_num[1])%cell_num[1]) * cell_num[0];
							
							projection_interaction_nsym();
						}
				}
		}

		if (phy_num_dim == 3)
		{
			for(cell_index[2]=0;cell_index[2]<cell_num[2];(cell_index[2])++)
				for(cell_index[1]=0;cell_index[1]<cell_num[1];(cell_index[1])++)
					for(cell_index[0]=0;cell_index[0]<cell_num[0];(cell_index[0])++)
					{
						cell_id = cell_index[0] + cell_index[1] * cell_num[0] + cell_index[2] * cell_num[1] * cell_num[0];
					
						for(neighbour_index[2]=cell_index[2]-1;neighbour_index[2]<=cell_index[2]+1;(neighbour_index[2])++)
							for(neighbour_index[1]=cell_index[1]-1;neighbour_index[1]<=cell_index[1]+1;(neighbour_index[1])++)	
								for(neighbour_index[0]=cell_index[0]-1;neighbour_index[0]<=cell_index[0]+1;(neighbour_index[0])++)
								{
									boundary_EDS();
								
									neighbour_id = ((neighbour_index[0] + cell_num[0])%cell_num[0]) + ((neighbour_index[1] + cell_num[1])%cell_num[1]) * cell_num[0] + ((neighbour_index[2] + cell_num[2])%cell_num[2])  * cell_num[1] * cell_num[0];

									if ((ctr_kernel_type == 3) || (ctr_kernel_type == 4))    
									    projection_interaction_remesh();

									else
									    projection_interaction_nsym();
									  
								}
					}
		}

		  /*
		  for(id=0;id<particles_num;id++)
		    for(ndim=0;ndim<phy_num_dim;ndim++)
		    {
			 particles_v_projection[id][ndim] = particles_v_projection[id][ndim] / particles_w[id];
			 particles_a_projection[id][ndim] = particles_a_projection[id][ndim] / particles_w[id];
		    }
		  */
	}
}

void ZS_Projection::projection_parameter_mb()
{
	//w = get_kernel(0);

	//for(id=0;id<particles_num;id++)
		//particles_rho_projection[id] = particles_rho_projection[id] + particles_m_projection[id] * w;

	if (ctr_list_type == 2)
	{
		if (phy_num_dim == 2)
		{
			for(cell_index[1]=0;cell_index[1]<cell_num[1];(cell_index[1])++)
				for(cell_index[0]=0;cell_index[0]<cell_num[0];(cell_index[0])++)
				{
					cell_id = cell_index[0] + cell_index[1] * cell_num[0];

					for(neighbour_index[1]=cell_index[1]-1;neighbour_index[1]<=cell_index[1]+1;(neighbour_index[1])++)
					{	
						if (neighbour_index[1] >= 0 && neighbour_index[1] < cell_num[1])
						{
							for(neighbour_index[0]=cell_index[0]-1;neighbour_index[0]<=cell_index[0]+1;(neighbour_index[0])++)
							{
								boundary_EDS();
								
								neighbour_id = ((neighbour_index[0] + cell_num[0])%cell_num[0])  + ((neighbour_index[1] + cell_num[1])%cell_num[1]) * cell_num[0];
								
								projection_interaction_nsym();
							}
						}
						else
						{
							boundary_lees_cell();

							for(neighbour_index[0]=le_min_index;neighbour_index[0]<=le_max_index;(neighbour_index[0])++)
							{
								if ((le_lower_cell_min == -999999) || (le_upper_cell_max == -999999))
									boundary_EDS();
								else
									boundary_NEDS();

								neighbour_id = ((neighbour_index[0] + cell_num[0])%cell_num[0])  + ((neighbour_index[1] + cell_num[1])%cell_num[1]) * cell_num[0];

								projection_interaction_nsym();
							}
						}
					}
				}
		}


		if (phy_num_dim == 3)
		{
			for(cell_index[2]=0;cell_index[2]<cell_num[2];(cell_index[2])++)
				for(cell_index[1]=0;cell_index[1]<cell_num[1];(cell_index[1])++)
					for(cell_index[0]=0;cell_index[0]<cell_num[0];(cell_index[0])++)
					{
						cell_id = cell_index[0] + cell_index[1] * cell_num[0] + cell_index[2] * cell_num[1] * cell_num[0];
					
						for(neighbour_index[2]=cell_index[2]-1;neighbour_index[2]<=cell_index[2]+1;(neighbour_index[2])++)
						{	
							if (neighbour_index[2] >= 0 && neighbour_index[2] < cell_num[2])
							{
								for(neighbour_index[1]=cell_index[1]-1;neighbour_index[1]<=cell_index[1]+1;(neighbour_index[1])++)	
									for(neighbour_index[0]=cell_index[0]-1;neighbour_index[0]<=cell_index[0]+1;(neighbour_index[0])++)
									{
										boundary_EDS();
								
										neighbour_id = ((neighbour_index[0] + cell_num[0])%cell_num[0]) + ((neighbour_index[1] + cell_num[1])%cell_num[1]) * cell_num[0] + ((neighbour_index[2] + cell_num[2])%cell_num[2])  * cell_num[1] * cell_num[0];

								
										projection_interaction_nsym();
									}
							}
							else
							{
								boundary_lees_cell();

								for(neighbour_index[1]=cell_index[1]-1;neighbour_index[1]<=cell_index[1]+1;(neighbour_index[1])++)
									for(neighbour_index[0]=le_min_index;neighbour_index[0]<=le_max_index;(neighbour_index[0])++)
									{
										if ((le_lower_cell_min == -999999) || (le_upper_cell_max == -999999))
											boundary_EDS();
										else
											boundary_NEDS();

										neighbour_id = ((neighbour_index[0] + cell_num[0])%cell_num[0])  + ((neighbour_index[1] + cell_num[1])%cell_num[1]) * cell_num[0];

										projection_interaction_nsym();
									}
							}
						}
					}
			}
	}
}


void ZS_Projection::projection_interaction_nsym()
{
	int ndim(0),num_i(0),num_j(0);
	double w(0.0),dr(0.0),r2_ij(0.0),r_ij(0.0);
	long i(0),j(0);
	vector<double> v_j;

	for(num_i=0;num_i<cell_list_projection[cell_id].size();num_i++)
	{
		i = cell_list_projection[cell_id][num_i];

		for(num_j=0;num_j<cell_list[neighbour_id].size();num_j++)
		{
			j = cell_list[neighbour_id][num_j];
		
			r2_ij = 0.0;

			for(ndim=0;ndim<phy_num_dim;ndim++)
			{
				dr = particles_r_projection[i][ndim] - (particles_r[j][ndim] + boundary_rshift[ndim]);
				r2_ij = r2_ij + dr * dr;
			}
			
			r_ij = sqrt(r2_ij);
			  //  cout << phy_cut_off << endl;
			
			if (r_ij < phy_cut_off)
			{
				w = get_kernel(r_ij);

				particles_w[i] = particles_w[i] + particles_m[j] / particles_rho[j] * w;

				if (ctr_only_velocity == 1)
					particles_rho_projection[i] = 0;
				else
					particles_rho_projection[i] = particles_rho_projection[i] + particles_m[j] * w;

				v_j = boundary_lees_velocity(j);

				for(ndim=0;ndim<phy_num_dim;ndim++)
				{
					particles_v_projection[i][ndim] = particles_v_projection[i][ndim] + particles_m[j] / particles_rho[j] * w * v_j[ndim];

					if (ctr_only_velocity == 1)
						particles_a_projection[i][ndim] = 0;
					else
						particles_a_projection[i][ndim] = particles_a_projection[i][ndim] + particles_m[j] / particles_rho[j] * w * particles_a[j][ndim];
				}	
			}
				
		}
	}
}

void ZS_Projection::projection_interaction_remesh()
{
        int ndim(0),num_i(0),num_j(0),r_flag(0);
	double w(1.0),w_temp(0.0);
	long i(0),j(0);
	vector<double> v_j;
	vector<double> dr(phy_num_dim, 0.0);

	for(num_i=0;num_i<cell_list_projection[cell_id].size();num_i++)
	{
		i = cell_list_projection[cell_id][num_i];

		for(num_j=0;num_j<cell_list[neighbour_id].size();num_j++)
		{
			j = cell_list[neighbour_id][num_j];
		
			w = 1.0;

			r_flag = 0;

			for(ndim=0;ndim<phy_num_dim;ndim++)
			{    
			    dr[ndim] = abs(particles_r_projection[i][ndim] - (particles_r[j][ndim] + boundary_rshift[ndim]));			    
			    
			    if (dr[ndim] > phy_cut_off) 
			        r_flag = 1;
			}

			if (r_flag == 0)
			{

			    for(ndim=0;ndim<phy_num_dim;ndim++)
			    {  
			        w_temp = get_kernel(dr[ndim]);

				w = w * w_temp;
			    }		

			    v_j = boundary_lees_velocity(j);

			    particles_rho_projection[i] = particles_rho_projection[i] + w * particles_rho[j];

			    for(ndim=0;ndim<phy_num_dim;ndim++)
			    {               			
					particles_v_projection[i][ndim] = particles_v_projection[i][ndim] + w * v_j[ndim];

					if (ctr_only_velocity == 1)
						particles_a_projection[i][ndim] = 0;
					else
						particles_a_projection[i][ndim] = particles_a_projection[i][ndim] + w * particles_a[j][ndim];
			    }	
			}
				
		}
	}
}


void ZS_Projection::projection_parameter_reset()
{
	long id(0);
	int ndim(0);

	for(id=0;id<particles_num;id++)
	{
		
		
		for(ndim=0;ndim<phy_num_dim;ndim++)
		{
			particles_v_projection[id][ndim] = 0.0;
			particles_a_projection[id][ndim] = 0.0;
		}
		
		particles_rho_projection[id] = 0.0;
		particles_w[id] = 0.0;
	}

	  //particles_r_projection.resize(particles_num);
}
