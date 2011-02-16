#pragma once
#include "ZS_Kernel.h"

class ZS_NeighbourList : public ZS_Kernel
{
	public:

		//vector< vector<long> >	neighbour_list;
		
		vector<long>			neighbour_count;
		vector<long>			neighbour_index;
		long					neighbour_id;
		
		// Cell List Method Parameter

		vector< vector<long> >	cell_list;
		vector< vector<long> >  cell_list_projection;
		vector<long>			cell_num;
		vector<long>			cell_index;
		vector<double>			cell_size;
		long					cell_id;
		long					ncells;

	public:
		
		ZS_NeighbourList();
		~ZS_NeighbourList();

		void cell_list_init();
		void cell_list_generate();
};
