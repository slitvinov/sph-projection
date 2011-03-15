///\file ZS_Boundary.h
#include "ZS_Boundary.h"

ZS_Boundary::ZS_Boundary()
{
}

ZS_Boundary::~ZS_Boundary()
{
}

void ZS_Boundary::boundary_init()
{
  boundary_upper_neighbour.resize(cell_num[0]+3, 0);
  boundary_lower_neighbour.resize(cell_num[0]+3, 0);
    boundary_rshift.resize(phy_num_dim, 0.0);
}

void ZS_Boundary::boundary_generate()
{
	if (ctr_test_case == 1)
		boundary_normal();

	if (ctr_test_case == 2)
		boundary_lees_edwards();

	if (ctr_test_case == 3)
		boundary_normal();

	if (ctr_test_case == 4)
		boundary_normal();

	if (ctr_test_case == 5)
		boundary_normal();

	if (ctr_test_case == 6)
		boundary_normal();
}

void ZS_Boundary::boundary_EDS()
{
	int ndim(0);

	for(ndim=0;ndim<phy_num_dim;ndim++)
	{
		if (neighbour_index[ndim] < 0)
			boundary_rshift[ndim] = - phy_domain_length[ndim];
		else if (neighbour_index[ndim] >= cell_num[ndim])
			boundary_rshift[ndim] = phy_domain_length[ndim];
		else
			boundary_rshift[ndim] = 0;
	}

}

void ZS_Boundary::boundary_NEDS()
{
	int ndim(0),ndim_vertical(0);

	ndim_vertical = phy_num_dim - 1;

	for(ndim=0;ndim<phy_num_dim;ndim++)
	{
		if (ndim == 0)
		{
			if (neighbour_index[ndim_vertical] < 0)
			{
				if (neighbour_index[0] < ncells - cell_num[0])
					boundary_rshift[0] = phy_domain_length[0] * (- fmod((phy_gamma_dot * phy_domain_length[ndim_vertical] * phy_t),phy_domain_length[0]) - 1);
				else if (neighbour_index[0] >= ncells)
					boundary_rshift[0] = phy_domain_length[0] * (- fmod((phy_gamma_dot * phy_domain_length[ndim_vertical] * phy_t),phy_domain_length[0]) + 1);
				else
					boundary_rshift[0] = phy_domain_length[0] * (- fmod((phy_gamma_dot * phy_domain_length[ndim_vertical] * phy_t),phy_domain_length[0]));
			}
			else if (neighbour_index[ndim_vertical] >= cell_num[ndim_vertical])
			{
				if (neighbour_index[0] < 0)
					boundary_rshift[0] = phy_domain_length[0] * (fmod((phy_gamma_dot * phy_domain_length[ndim_vertical] * phy_t),phy_domain_length[0]) - 1);				
				else if (neighbour_index[0] >= cell_num[0])
					boundary_rshift[0] = phy_domain_length[0] * (fmod((phy_gamma_dot * phy_domain_length[ndim_vertical] * phy_t),phy_domain_length[0]) + 1);
				else
					boundary_rshift[0] = phy_domain_length[0] * (fmod((phy_gamma_dot * phy_domain_length[ndim_vertical] * phy_t),phy_domain_length[0]));
			}
			else
			{
				if (neighbour_index[ndim] < 0)
					boundary_rshift[ndim] = - phy_domain_length[ndim];
				else if (neighbour_index[ndim] >= cell_num[ndim])
					boundary_rshift[ndim] = phy_domain_length[ndim];
				else
					boundary_rshift[ndim] = 0;
			}	
		}

		else
		{
			if (neighbour_index[ndim] < 0)
				boundary_rshift[ndim] = - phy_domain_length[ndim];
			else if (neighbour_index[ndim] >= cell_num[ndim])
				boundary_rshift[ndim] = phy_domain_length[ndim];
			else
				boundary_rshift[ndim] = 0;
		}
	}
}


void ZS_Boundary::boundary_lees_edwards()
{
	double x_upper(0.0),x_lower(0.0);
	long i(0);
	int ndim_vertical(0);

	ndim_vertical = phy_num_dim - 1;

	x_upper = phy_domain_length[0] * (1 - fmod((phy_gamma_dot * phy_domain_length[ndim_vertical] * phy_t),phy_domain_length[0]));
	x_lower = phy_domain_length[0] * (1 - fmod((phy_gamma_dot * phy_domain_length[ndim_vertical] * phy_t),phy_domain_length[0]));

	if (x_upper == phy_domain_length[0])	
		le_upper_cell_max = -999999;
	else
	{
		le_upper_cell_max = ceil(x_upper / cell_size[0]);

		if ((le_upper_cell_max - (cell_num[0] + 2)) < (-cell_num[0]))
			le_upper_cell_max = le_upper_cell_max + cell_num[0];

		for(i=0;i<cell_num[0]+3;i++)
			boundary_upper_neighbour[i] = le_upper_cell_max - (cell_num[0] + 2 - i);
	}


	if (x_lower == phy_domain_length[0])
		le_lower_cell_min = -999999;
	else
	{
		le_lower_cell_min = (ncells - 1) - ceil(x_lower / cell_size[0]);

		for(i=0;i<cell_num[0]+3;i++)
			boundary_lower_neighbour[i] = le_lower_cell_min + i;
	}		

}

void ZS_Boundary::boundary_normal()
{
	le_upper_cell_max = -999999;
	le_lower_cell_min = -999999;
}

void ZS_Boundary::boundary_lees_cell()
{
	int ndim_vertical(0);

	ndim_vertical = phy_num_dim - 1;

	if (neighbour_index[ndim_vertical] < 0)
	{
		if (le_lower_cell_min == -999999)
		{
			le_min_index = cell_index[0]-1;
			le_max_index = cell_index[0]+1;	
		}
		else
		{
			le_min_index = boundary_lower_neighbour[cell_index[0]];
			le_max_index = boundary_lower_neighbour[cell_index[0]+3];
		}
	}
	else
	{
		if ( le_upper_cell_max == -999999)
		{	
			le_min_index = cell_index[0]-1;
			le_max_index = cell_index[0]+1;
		}
		else
		{
			le_min_index = boundary_upper_neighbour[cell_index[0]];
			le_max_index = boundary_upper_neighbour[cell_index[0]+3];
		}
	}
}

vector<double> ZS_Boundary::boundary_lees_velocity(long j)
{
	int ndim(0),ndim_vertical(0);
	vector<double> v_j(phy_num_dim,0.0);

	ndim_vertical = phy_num_dim - 1;

	if (neighbour_index[ndim_vertical] < 0)
	{
		v_j[0] = (particles_v[j][0] - phy_gamma_dot * phy_domain_length[ndim_vertical]);
		for(ndim=1;ndim<phy_num_dim;ndim++)
			v_j[ndim] = particles_v[j][ndim];
	}
	else if (neighbour_index[ndim_vertical] >= cell_num[ndim_vertical])
	{
		v_j[0] = (particles_v[j][0] + phy_gamma_dot * phy_domain_length[ndim_vertical]);
		for(ndim=1;ndim<phy_num_dim;ndim++)
			v_j[ndim] = particles_v[j][ndim];
	}
	else
	{
		for(ndim=0;ndim<phy_num_dim;ndim++)
			v_j[ndim] = particles_v[j][ndim];
	}

	return v_j;
}
