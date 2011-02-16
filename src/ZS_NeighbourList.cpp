#include "ZS_NeighbourList.h"

ZS_NeighbourList::ZS_NeighbourList()
{		
}

ZS_NeighbourList::~ZS_NeighbourList()
{
}

void ZS_NeighbourList::cell_list_init()
{
	long i(0),index(0);
	double A_err(10e-15);
	int ndim(0);

	ncells = 1;

	for(i=0;i < phy_num_dim;i++)
	{
		cell_num.push_back(floor(phy_domain_length[i]/ phy_cut_off));
		cell_size.push_back(phy_domain_length[i] / cell_num[i]);
		ncells = ncells * cell_num[i];
	}

	  // cout << ncells << endl;
	  //cout << phy_cut_off << endl;
	cell_list.resize(ncells);
	cell_list_projection.resize(ncells);
	cell_index.resize(phy_num_dim,0);
	neighbour_index.resize(phy_num_dim,0);

	for(i=0;i<particles_num;i++)
	{
		//if (i == particles_num - 2) ;
		for(ndim=0;ndim<phy_num_dim;ndim++)
			cell_index[ndim] = floor(particles_r_projection[i][ndim] / (cell_size[ndim]+A_err));
		
		if (phy_num_dim == 2)
			cell_id = cell_index[0] + cell_index[1] * cell_num[0];

		if (phy_num_dim == 3)
			cell_id = cell_index[0] + cell_index[1] * cell_num[0] + cell_index[2] * cell_num[0] * cell_num[1];

		  //if ((cell_id<1000) && (cell_id >= 0))
			cell_list_projection[cell_id].push_back(i);
		  //else 
			  //{	index = i;
		    //cout << "partilce " << i << endl;
			  //}
	}
}

void ZS_NeighbourList::cell_list_generate()
{
	long i(0),index(0);
	int ndim(0);

	for(i=0;i<ncells;i++)
		cell_list[i].clear();

	for(i=0;i<particles_num;i++)
	{
		for(ndim=0;ndim<phy_num_dim;ndim++)
			cell_index[ndim] = floor(particles_r[i][ndim] / cell_size[ndim]);
		
		if (phy_num_dim == 2)
			cell_id = cell_index[0] + cell_index[1] * cell_num[0];

		if (phy_num_dim == 3)
			cell_id = cell_index[0] + cell_index[1] * cell_num[0] + cell_index[2] * cell_num[0] * cell_num[1];

		//if ((cell_id<1000) && (cell_id>= 0))
			cell_list[cell_id].push_back(i);
		//else
		//	index = i;
	}
}
