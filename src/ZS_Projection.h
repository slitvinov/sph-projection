#pragma once
#include "ZS_Boundary.h"

class ZS_Projection : public ZS_Boundary
{
	public:
	
		ZS_Projection();
		~ZS_Projection();

		void projection_init();
		void projection_calculate();
		void projection_position();
		void projection_parameter();
		void projection_parameter_pb();
		void projection_parameter_mb();
		void projection_interaction_nsym();
		void projection_interaction_remesh();
		void projection_parameter_reset();
};
