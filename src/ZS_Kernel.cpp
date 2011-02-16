#include "ZS_Kernel.h"

ZS_Kernel::ZS_Kernel()
{
}

ZS_Kernel::~ZS_Kernel()
{
}

double ZS_Kernel::get_kernel(double r)
{
	double s(0.0);
	double factor(0.0);
	  double w(0.0);

	// Piecewise Quintic Spline	(Morris, 1996)
	if (ctr_kernel_type == 1)
	{
		s = 3 * r / (phy_cut_off);

		if (phy_num_dim == 3)    
			factor = 27 / (120*pi*pow(phy_cut_off,3)) ;
				
		if (phy_num_dim == 2)    
			factor = 63 / (478*pi*pow(phy_cut_off,2)) ;

		if (s>=0 && s<1) 
			w = factor * ( pow((3-s),5) - 6*pow((2-s),5) + 15*pow((1-s),5) ) ;
		if (s>=1 && s<2)
			w = factor * ( pow((3-s),5) - 6*pow((2-s),5) ) ;
		if (s>=2 && s<3)
			w = factor *   pow((3-s),5) ;
		if (s<0 || s>=3)
			w = 0 ;
	}

	// Quartic Spline	(Lucy, 1977)
	if (ctr_kernel_type == 2)
	{
		s = r / phy_cut_off;

		if (phy_num_dim == 3)
			factor = 105 / (16*pi*pow(phy_cut_off,3));

		if (phy_num_dim == 2)
			factor = 5 / (pi*pow(phy_cut_off,2));

		if (s<=1)
			w = factor * (1+3*s) * pow((1-s),3);
		else
			w = 0;
	}

	if (ctr_kernel_type == 3)
	{	       
	        s = 2 * r / phy_cut_off;

		if (s>=0 && s<1)
			w = 1 - 5 * s * s / 2 + 3 * s * s * s / 2;
		if (s>=1 && s<2)
		        w = (1-s)*(2-s)*(2-s) / 2;
		if (s>=2)
			w = 0;
	}

	if (ctr_kernel_type == 4)
	{	       
	        s = 3 * r / phy_cut_off;

		if (s>=0 && s<1)
		        w = (1-s) * (25*pow(s,4) - 38*pow(s,3) - 3*s*s + 12*s + 12) / 12;
		if (s>=1 && s<2)
		        w = (s-1) * (s-2) * (25*pow(s,3) - 114*s*s + 153*s - 48) / 24;
		if (s>=2 && s<3)
		        w = (2-s) * pow((s-3),3) * (5*s - 8) / 24; 
		if (s>=3)
			w = 0;
	}

	return w;
}

double ZS_Kernel::get_kernel_derivity(double r)
{
	double s(0.0);
	double factor(0.0);
	double gradw(0.0);

	if (ctr_kernel_type == 1)
	{
		s = 3 * r / (phy_cut_off);

		if (phy_num_dim == 3)    
			factor = 27 / (120*pi*pow(phy_cut_off,3)) ;
				
		if (phy_num_dim == 2)    
			factor = 63 / (478*pi*pow(phy_cut_off,2)) ;

		if (s>=0 && s<1) 
			gradw = factor * (-15) / phy_cut_off * ( pow((3-s),4) - 6*pow((2-s),4) + 15*pow((1-s),4) ) ;
		if (s>=1 && s<2)
			gradw = factor * (-15) / phy_cut_off * ( pow((3-s),4) - 6*pow((2-s),4) ) ;
		if (s>=2 && s<3)
			gradw = factor * (-15) / phy_cut_off * ( pow((3-s),4) ) ;
		if (s<0 || s>=3)
			gradw = 0 ;
	}

	if (ctr_kernel_type == 2)
	{
		s = r / phy_cut_off;

		if (phy_num_dim == 3)
			factor = 105 / (16*pi*pow(phy_cut_off,3));

		if (phy_num_dim == 2)
			factor = 5 / (pi*pow(phy_cut_off,2));

		if (s<=1)
			gradw = factor / phy_cut_off * ( 3 * pow((1-s),3) - 3 * (1+3*s) * pow((1-s),2) ); 
		else
			gradw = 0;
	}


	return gradw;
}
