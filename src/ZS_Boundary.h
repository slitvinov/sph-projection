#ifndef ZS_BOUNDARY_H
#define ZS_BOUNDARY_H
#include "ZS_NeighbourList.h"

class ZS_Boundary : public ZS_NeighbourList
{
	public:

			long		le_upper_cell_max;
			long		le_lower_cell_min;
			long		le_max_index;
			long		le_min_index;

		vector<long>	boundary_upper_neighbour;
		vector<long>	boundary_lower_neighbour;
		vector<double>	boundary_rshift;

	public:

		ZS_Boundary();
		~ZS_Boundary();

		void boundary_init();
		void boundary_generate();

		void boundary_EPBC(long id,int ndim);
		void boundary_NEPBC(long id,int ndim);

		void boundary_EDS();
		void boundary_NEDS();

		void boundary_normal();

		void boundary_lees_edwards();
		void boundary_lees_cell();
		vector<double> boundary_lees_velocity(long j);
};

#endif //ZS_BOUNDARY_H
