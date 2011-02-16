#include "ZS_Physics.h"

ZS_Physics::ZS_Physics()
{
}

ZS_Physics::~ZS_Physics()
{
}

void ZS_Physics::phy_init()
{
  phy_read();

       
  phy_domain_length.resize(phy_num_dim,0.0);

  //cout << "yes!!" << endl;
  get_domain_length();
  get_cut_off();
  get_time_step();

	 
  phy_momentum.resize(phy_num_dim,0.0);
  phy_kinetic_energy = 0.0;
  phy_random_seed = -time(NULL);
  srand(time(NULL));
  phy_zeta_max = 0.0;
  phy_zeta_min = 0.0;
	
}

void ZS_Physics::phy_read()
{
  string cbuf,carg,cvalue,trimmstr(" "),headstr("#"),equalstr("=");
  signed int loc;
  double value;
  vector<double> *vvalue =  new vector<double>();

  ifstream inputFile(ctr_physics_config_name.c_str());

  if( !inputFile ) 
  {
    cerr << "Can't open physics file" << endl;    
    return;
  }

  while( !inputFile.eof() ) 
  {
    getline(inputFile,cbuf);

    if (cbuf.size() <= 1 ) continue;
    else
    {		
      if (cbuf.find_first_of(headstr,0) != string::npos) continue;
      else 
      {	
        loc = cbuf.find(equalstr, 0);	
        if( loc == string::npos ) cout << "Error in physics file! " << endl;	
        carg=cbuf.substr(0,loc-1);	
        cvalue=cbuf.substr(loc+1);	

        // for linux version
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

        value = atof(cvalue.c_str());

        if (carg.compare("num_dim") == 0)
          phy_num_dim = int(value);					
        if (carg.compare("min_length") == 0)
        {
          phy_get_vector(cvalue,vvalue);
          phy_min_length.assign(vvalue->begin(),vvalue->end());
          vvalue->clear();
        }
        if (carg.compare("max_length") == 0)
        {
          phy_get_vector(cvalue,vvalue);
          phy_max_length.assign(vvalue->begin(),vvalue->end());
          vvalue->clear();
        }
        if (carg.compare("num_particles") == 0)
        {
          phy_get_vector(cvalue,vvalue);
          phy_num_particles.assign(vvalue->begin(),vvalue->end());
          vvalue->clear();
        }
        if (carg.compare("num_macroparticles") == 0)
        {
          phy_get_vector(cvalue,vvalue);
          phy_num_macroparticles.assign(vvalue->begin(),vvalue->end());
          vvalue->clear();
        }
        if (carg.compare("overlap") == 0)
          phy_overlap = value;
        if (carg.compare("dt") == 0)
          phy_dt = value;
        if (carg.compare("tau_dt") == 0)
          phy_tau_dt = value;
        if (carg.compare("tau_macro") == 0)
          phy_tau_macro = value;
        if (carg.compare("time_start") == 0)
          phy_time_start = value;
        if (carg.compare("time_end") == 0)
          phy_time_end= value;
        if (carg.compare("step_start") == 0)
          phy_step_start = long(value);
        if (carg.compare("step_end") == 0)
          phy_step_end = long(value);
        if (carg.compare("msd_step") == 0)
          phy_msd_step = long(value);
        if (carg.compare("v0") == 0)
          phy_v0 = value;
        if (carg.compare("v0_macro") == 0)
          phy_v0_macro = value;
        if (carg.compare("rho") == 0)
          phy_rho = value;
        if (carg.compare("eta") == 0)
          phy_eta = value;
        if (carg.compare("ksai") == 0)
          phy_ksai = value;
        if (carg.compare("temp") == 0)
          phy_temp = value;
        if (carg.compare("k0") == 0)
          phy_k0 = value;
        if (carg.compare("c") == 0)
          phy_c = value;
        if (carg.compare("rho_0") == 0)
          phy_rho_0 = value;
        if (carg.compare("rho_0_macro") == 0)
          phy_rho_0_macro = value;
        if (carg.compare("b") == 0)
          phy_b = value;
        if (carg.compare("gamma") == 0)
          phy_gamma = value;
        if (carg.compare("gamma_dot") == 0)
          phy_gamma_dot = value;
        if (carg.compare("Fc") == 0)
          phy_Fc = value;
        if (carg.compare("Fr") == 0)
          phy_Fr = value;
        if (carg.compare("epsilon_0") == 0)
          phy_epsilon_0 = value;
        if (carg.compare("Ma") == 0)
          phy_Ma = value;
        if (carg.compare("stochastic_A") == 0)
          phy_A = value;
        if (carg.compare("body_force") == 0)
        {
          phy_get_vector(cvalue,vvalue);
          phy_body_force.assign(vvalue->begin(),vvalue->end());
          vvalue->clear();
        }
        if (carg.compare("offset") == 0)
        {
          phy_get_vector(cvalue,vvalue);
          phy_offset.assign(vvalue->begin(),vvalue->end());
          vvalue->clear();
        }
        if (carg.compare("bcdef") == 0)
        {
          phy_get_vector(cvalue,vvalue);
          phy_bcdef.assign(vvalue->begin(),vvalue->end());
          vvalue->clear();
        }
      }
    }

  }

  inputFile.close();

  delete(vvalue);
}


void ZS_Physics::phy_print(FILE *fid)
{
  cout << "\n" << endl;
  fprintf(fid,"\n");
  cout << "#################################################" << endl;
  fprintf(fid,"#################################################\n");
  cout << "##                    ZS_SPH                   ##" << endl;
  fprintf(fid,"##                    ZS_SPH                   ##\n");
  cout << "##                 Version : 1.1               ##" << endl;
  fprintf(fid,"##                 Version : 1.1               ##\n");
  cout << "##                                             ##" << endl;
  fprintf(fid,"##                                             ##\n");
  cout << "##                 Autor : Yilei Shi           ##" << endl;
  fprintf(fid,"##                 Autor : Yilei Shi           ##\n");
  cout << "##                   2010.03.10                ##" << endl;
  fprintf(fid,"##                   2010.03.10                ##\n");
  cout << "#################################################" << endl;
  fprintf(fid,"#################################################\n");
  cout << "\n" << endl;
  fprintf(fid,"\n");
  cout << "###########     Physical Properties   ###########" << endl;
  fprintf(fid,"###########     Physical Properties   ###########\n");
  cout << "\n" << endl;
  fprintf(fid,"\n");

  if (ctr_test_case == 1)
    cout << " Test Case " << ctr_test_case << "                  :  Euler Flow" << endl;
  if (ctr_test_case == 2)
    cout << " Test Case " << ctr_test_case << "                  :  Shear Flow" << endl;
  if (ctr_test_case == 3)
    cout << " Test Case " << ctr_test_case << "                  :  Taylor-Green Vortex" << endl;
  if (ctr_test_case == 4)
    cout << " Test Case " << ctr_test_case << "                  :  Forced Turbulence" << endl;
  if (ctr_test_case == 5)
    cout << " Test Case " << ctr_test_case << "                  :  Stochastic Model (Two sets of particles)" << endl;
  if (ctr_test_case == 6)
    cout << " Test Case " << ctr_test_case << "                  :  Stochastic Model (One sets of particles)" << endl;
  cout << " Dimension                    :  D     =    " << phy_num_dim << endl;
  fprintf(fid," Dimension                    :  D     =     %d\n", phy_num_dim);
  if (phy_num_dim == 2)
  {
    cout << " Number of micro particles    :  N     =    " << phy_num_particles[0] << "," << phy_num_particles[1] << endl;
    fprintf(fid," Number of micro particles    :  N     =     %f,%f\n", phy_num_particles[0], phy_num_particles[1]);

    if (ctr_test_case == 5)
    {
      cout << " Number of macro particles    :  N     =    " << phy_num_macroparticles[0] << "," << phy_num_macroparticles[1] << endl;
      fprintf(fid," Number of macro particles    :  N     =     %f,%f\n", phy_num_macroparticles[0], phy_num_macroparticles[1]);
    }
    cout << " Length of domains            :  L     =    " << phy_domain_length[0] << "," << phy_domain_length[1] << endl;
    fprintf(fid," Length of domains            :  L     =     %f,%f\n", phy_domain_length[0], phy_domain_length[1]);
  }
  if (phy_num_dim == 3)
  {
    cout << " Number of micro particles    :  N     =    " << phy_num_particles[0] << "," << phy_num_particles[1] << "," << phy_num_particles[2] << endl;
    fprintf(fid," Number of micro particles    :  N     =     %f,%f,%f\n", phy_num_particles[0], phy_num_particles[1], phy_num_particles[2]);

    if (ctr_test_case == 5)
    {
      cout << " Number of macro particles    :  N     =    " << phy_num_macroparticles[0] << "," << phy_num_macroparticles[1] << "," << phy_num_macroparticles[2] << endl;
      fprintf(fid," Number of macro particles    :  N     =     %f,%f,%f\n", phy_num_macroparticles[0], phy_num_macroparticles[1], phy_num_macroparticles[2]);
    }
    cout << " Length of domains            :  L     =    " << phy_domain_length[0] << "," << phy_domain_length[1] << "," << phy_domain_length[2] << endl;
    fprintf(fid," Length of domains            :  L     =     %f,%f,%f\n", phy_domain_length[0], phy_domain_length[1], phy_domain_length[2]);
  }

  cout << " Viscosity                    :  nu    =    " << phy_eta << endl;
  fprintf(fid," Viscosity			     :  nu    =     %f\n", phy_eta);
  cout << " Overlap                      :  k     =    " << phy_cut_off * phy_num_particles[0] / phy_domain_length[0] << endl;
  fprintf(fid," Overlap                      :  k     =     %f\n", phy_cut_off * phy_num_particles[0] / phy_domain_length[0]);
  cout << " Smoothing length(micro)      :  h     =    " << phy_cut_off << endl;
  fprintf(fid," Smoothing length(micro)      :  h     =     %f\n", phy_cut_off);


  if (ctr_test_case == 5)
  {
    cout << " Smoothing length(macro)      :  h     =    " << phy_cut_off_macro << endl;
    fprintf(fid," Smoothing length(macro)      :  h     =     %f\n", phy_cut_off_macro);
  }
  cout << " Time step                    :  dt    =    " << phy_dt << endl;
  fprintf(fid," Time step                    :  dt    =     %f\n", phy_dt);
  cout << " Time step(CFL)               :  dt    =    " << phy_dt_CFL << endl;
  fprintf(fid," Time step(CFL)               :  dt    =     %f\n", phy_dt_CFL);
  cout << " Speed of sound               :  c     =    " << phy_c << endl;
  fprintf(fid," Speed of sound               :  c     =     %f\n", phy_c);

  if (ctr_test_case == 2)
    cout << " Tensor calculation           :   " << ctr_tensor_calculate << endl;

  if (ctr_test_case == 4)
  {
    cout << " Rotational force             :  Fr    =    " << phy_Fr << endl;
    fprintf(fid," Rotational force             :  Fr    =     %f\n", phy_Fr);
    cout << " Compressible force           :  Fc    =    " << phy_Fc << endl;
    fprintf(fid," Compressible force           :  Fc    =     %f\n", phy_Fc);
    cout << " Dissipation rate             :  e     =    " << phy_epsilon << endl;
    fprintf(fid," Dissipation rate			   :  e     =     %f\n", phy_epsilon);
    cout << " Kolmogorov scale             :  eta   =    " << phy_kolm_eta << endl;
    fprintf(fid," Kolmogorov scale             :  eta   =     %f\n", phy_kolm_eta);
    cout << " Kolmogorov time              :  tau   =    " << phy_kolm_tau << endl;
    fprintf(fid," Kolmogorov time              :  tau   =     %f\n", phy_kolm_tau);	
    cout << " Before external force        :   " << ctr_bfext << endl;
  }

  if (ctr_test_case == 5)
  {
    cout << " Relaxation time              :  tau   =    " << phy_tau_macro << endl;
    fprintf(fid," Relaxation time              :  tau   =     %f\n", phy_tau_macro);
  }

  fclose(fid);

  cout << " Press any Key to continue... " << endl;
  //cbuf = getchar();
}


void ZS_Physics::phy_get_vector(string cvector,vector<double> *vector_value)		
{
  string gcbuf,gcvalue,gctemp,trimstr(",");
  int i(0);
  signed int loc;

  gcbuf = cvector;
  loc = gcbuf.find(trimstr, 0);

  while (loc != string::npos)
  {
    gcvalue=gcbuf.substr(0,loc);

    //* for linux version
    gcvalue.erase(remove(gcvalue.begin(), gcvalue.end(), ' '), gcvalue.end());		
    gcvalue.erase(remove(gcvalue.begin(), gcvalue.end(), '\t'), gcvalue.end());		
    gcvalue.erase(remove(gcvalue.begin(), gcvalue.end(), '\n'), gcvalue.end());
		

    //This line is for windows version
    //gcvalue.erase(remove_if(gcvalue.begin(), gcvalue.end(), isspace), gcvalue.end());

    vector_value->push_back(atof(gcvalue.c_str()));
    gctemp = gcbuf.substr(loc+1);
    gcbuf.assign(gctemp);
    loc = gcbuf.find(trimstr, 0);

    i++;
  }

  if (i == 0) cout << "Error in control file! " << endl;

  // for linux version
  gcbuf.erase(remove(gcbuf.begin(), gcbuf.end(), ' '), gcbuf.end());
  gcbuf.erase(remove(gcbuf.begin(), gcbuf.end(), '\t'), gcbuf.end());
  gcbuf.erase(remove(gcbuf.begin(), gcbuf.end(), '\n'), gcbuf.end());
	

  //This line is for windows version
  //gcbuf.erase(remove_if(gcbuf.begin(), gcbuf.end(), isspace), gcbuf.end());

  vector_value->push_back(atof(gcvalue.c_str()));
}


void ZS_Physics::get_domain_length()
{
  int ndim(0);

  for(ndim=0;ndim<phy_num_dim;ndim++)
    phy_domain_length[ndim] = phy_max_length[ndim] - phy_min_length[ndim];
}

void ZS_Physics::get_time_step()
{
  double dt(0.0),v(0.0);

  if (phy_c >= phy_v0)	v = phy_c;
  else	v = phy_v0;

  dt = 0.1 * 0.25 * phy_cut_off / v ;


  phy_dt_CFL = dt;	

  if (phy_eta != 0)
    dt =  0.1 * 0.125 * pow(phy_cut_off,2) / (phy_eta / phy_rho_0);
	
  //cout << "dt_eta :" << dt << endl;    
  if (dt < phy_dt_CFL) phy_dt_CFL = dt;	

  if (ctr_test_case == 4)
  {
    if (phy_epsilon_0 != -1)
      phy_epsilon = phy_epsilon_0;
    else
      phy_epsilon = phy_rho_0 * (phy_Fr + phy_Fc) ;

    phy_kolm_eta = double(pow(phy_eta*phy_eta*phy_eta / phy_epsilon, 0.25));
    phy_kolm_tau = sqrt(phy_eta / phy_epsilon);
    phy_dt = phy_tau_dt * phy_kolm_tau;

    if (ctr_dt_auto == 1)
    {
      if (phy_dt > phy_dt_CFL)
      {
        phy_dt = phy_dt_CFL;
        cout << "\n" << endl;
        cout << "Caution !!! : CFL condition is not satisfied! Chang small time step!" << endl;
      }
    }
    else
      cout << "Caution !!! : CFL condition is not satisfied! Chang small time step!" << endl;
  }

  else
  {
    if (ctr_dt_auto == 1)
      phy_dt = phy_dt_CFL;
    else
      cout << "Caution !!! : CFL condition is not satisfied! manuell time step!" << endl;
  }
  //        	    cout << "phy_dt :" << phy_dt << endl;  
}

void ZS_Physics::get_cut_off()
{
  if (ctr_kernel_type == 3)
    phy_overlap = 2;
  if (ctr_kernel_type == 4)
    phy_overlap = 3;
	
  phy_cut_off = phy_overlap * (phy_max_length[0] - phy_min_length[0]) / phy_num_particles[0];
  phy_cut_off_macro = phy_overlap * (phy_max_length[0] - phy_min_length[0]) / phy_num_macroparticles[0];
}
