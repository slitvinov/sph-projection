#include "ZS_Control.h"

ZS_Control::ZS_Control()
{
  ctr_io_name = "ctrl.sz";
}

ZS_Control::~ZS_Control()
{
}

void ZS_Control::ctrl_read()
{
  string cbuf,carg,cvalue;
  const string trimmstr(" ");
  const string headstr("#");
  const string equalstr("=");
  signed int loc;
  int value;

  ifstream inputFile(ctr_io_name.c_str());

	 
  if( !inputFile ) 
  {
    cerr << "Can't open control file" << endl;
    return;
  }

  while( !inputFile.eof() ) 
  {
    getline(inputFile,cbuf);

    if (cbuf.size() <= 1) ;
    else
    {		
      if (cbuf.find_first_of(headstr,0) != string::npos) ;
      else 
      {	
        loc = cbuf.find(equalstr, 0);	
        if( loc == string::npos ) cout << "Error in control file! " << endl;
        carg=cbuf.substr(0,loc-1);	
        cvalue=cbuf.substr(loc+1);

        //  for linux version
        carg.erase(remove(carg.begin(), carg.end(), ' '), carg.end());		
        carg.erase(remove(carg.begin(), carg.end(), '\t'), carg.end());		
        carg.erase(remove(carg.begin(), carg.end(), '\n'), carg.end());			
        cvalue.erase(remove(cvalue.begin(), cvalue.end(), ' '), cvalue.end());
        cvalue.erase(remove(cvalue.begin(), cvalue.end(), '\t'), cvalue.end());
        cvalue.erase(remove(cvalue.begin(), cvalue.end(), '\n'), cvalue.end());
        cvalue.erase(cvalue.length()-1, 1);
				
				
        // These two lines are for windows version 
        //carg.erase(remove_if(carg.begin(), carg.end(), isspace), carg.end());		
        //cvalue.erase(remove_if(cvalue.begin(), cvalue.end(), isspace), cvalue.end());
				
        value = atoi(cvalue.c_str());

        if (carg.compare("PHYSICS_CONFIG_FILE") == 0)
          ctr_physics_config_name = cvalue;
        if (carg.compare("READ_INIT_PARTICLES") == 0)
          ctr_read_init_particles = value;
        if (carg.compare("READ_INIT_MACRO_PARTICLES") == 0)
          ctr_read_init_macroparticles = value;
        if (carg.compare("INPUT_NAME") == 0)
          ctr_input_name = cvalue;
        if (carg.compare("INPUT_NAME_MACRO") == 0)
          ctr_input_name_macro = cvalue;
        if (carg.compare("INPUT_FORMAT") == 0)
          ctr_input_fmt = cvalue;
        if (carg.compare("OUTPUT_NAME") == 0)
          ctr_output_name = cvalue;
        if (carg.compare("OUTPUT_NAME_MACRO") == 0)
          ctr_output_name_macro = cvalue;
        if (carg.compare("OUTPUT_FORMAT") == 0)
          ctr_output_fmt = cvalue;
        if (carg.compare("OUTPUT_FREQUENCY") == 0)
          ctr_output_freq = value;
        if (carg.compare("STATISTIC_NAME") == 0)
          ctr_statistic_name = cvalue;
        if (carg.compare("STATISTIC_FORMAT") == 0)
          ctr_statistic_fmt = cvalue;
        if (carg.compare("STATISTIC_FREQUENCY") == 0)
          ctr_statistic_freq = value;
        if (carg.compare("WRITE_RESTART") == 0)
          ctr_write_restart = value;
        if (carg.compare("RESTART_NAME") == 0)
          ctr_restart_name = cvalue;
        if (carg.compare("RESTART_FORMAT") == 0)
          ctr_restart_fmt = cvalue;
        if (carg.compare("RESTART_FREQUENCY") == 0)
          ctr_restart_freq = value;
        if (carg.compare("DEBUG_FLAG") == 0)
          ctr_debug_flag = value;
        if (carg.compare("RHS_PRESSURE_FORM") == 0)
          ctr_rhs_pressure_form = value;
        if (carg.compare("RHS_VISCOSITY_FORM") == 0)
          ctr_rhs_viscosity_form = value;
        if (carg.compare("RHS_EXTERNAL_FORM") == 0)
          ctr_rhs_external_form = value;
        if (carg.compare("RHS_EXTERNAL_TYPE") == 0)
          ctr_rhs_external_type = value;
        if (carg.compare("STATE_EQUATION_TYPE") == 0)
          ctr_state_equation_type = value;
        if (carg.compare("STATE_EQUATION_FORM") == 0)
          ctr_state_equation_form = value;
        if (carg.compare("KERNEL_TYPE") == 0)
          ctr_kernel_type = value;
        if (carg.compare("NEIGHBOUR_LIST_TYPE") == 0)
          ctr_list_type = value;
        if (carg.compare("INTEGRATE_TYPE") == 0)
          ctr_integrate_type = value;
        if (carg.compare("SYMMETRY") == 0)
          ctr_symmetry = value;
        if (carg.compare("TEST_CASE") == 0)
          ctr_test_case = value;
        if (carg.compare("BEFORE_EXTERNAL") == 0)
          ctr_bfext = value;
        if (carg.compare("RANDOM_POSITION") == 0)
          ctr_random_position = value;
        if (carg.compare("RANDOM_VELOCITY") == 0)
          ctr_random_velocity = value;
        if (carg.compare("TENSOR_CALCULATE") == 0)
          ctr_tensor_calculate = value;
        if (carg.compare("SF_CALCULATE") == 0)
          ctr_sf_calculate = value;
        if (carg.compare("Auto_Time_Step") == 0)
          ctr_dt_auto = value;
        if (carg.compare("PROJECTION_TYPE") == 0)
          ctr_projection_type = value;
        if (carg.compare("ONLY_VELOCITY") == 0)
          ctr_only_velocity = value;
      }
    }

  }

  inputFile.close();
}
