#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <algorithm>
#include <cmath>
#include <ctime>
//#include <omp.h>

#define pi 3.1415926535897932384626

using namespace std;

class ZS_Control
{
	public:

			int			ctr_debug_flag;
			int			ctr_rhs_pressure_form;
			int			ctr_rhs_viscosity_form;
			int			ctr_rhs_external_form;
			int                     ctr_rhs_external_type;
                        int			ctr_state_equation_type;
			int			ctr_state_equation_form;
                        int			ctr_kernel_type;
			int			ctr_list_type;
                        int			ctr_integrate_type;
                        int			ctr_symmetry;
                        int			ctr_read_init_particles;
			int                     ctr_read_init_macroparticles;
                        int			ctr_write_restart;
			int		        ctr_output_freq;
			int			ctr_restart_freq;
			int			ctr_statistic_freq;
			int			ctr_test_case;
		        int                     ctr_bfext; 
			int			ctr_random_position;
			int			ctr_random_velocity;
			int			ctr_tensor_calculate;
                        int			ctr_sf_calculate;
			int                     ctr_dt_auto;
			int                     ctr_projection_type;
			int                     ctr_only_velocity;
			
			string		ctr_input_name;
			string		ctr_input_name_macro;
			string		ctr_input_fmt;
			string		ctr_output_name;
			string		ctr_output_name_macro;
			string		ctr_output_fmt;
			string		ctr_statistic_name;
			string		ctr_statistic_fmt;
			string		ctr_restart_name;
			string		ctr_restart_fmt;
			string		ctr_physics_config_name;
			string		ctr_io_name;

	public:
	
		ZS_Control();
		~ZS_Control();

		void ctrl_read();
};
