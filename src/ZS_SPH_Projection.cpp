#include "ZS_Projection.h"

int main(int argc,char** argv)
{
	char c;
	time_t time_start,time_end;
	double i(0);

	ZS_Projection sph_projection;

	time(&time_start);

	sph_projection.projection_init();

	  //cout << "right!!!" << endl;

	sph_projection.projection_calculate();

	time(&time_end);

	cout << "Caculation Time : " << difftime(time_end,time_start) <<  endl;

            //	c=getchar();

	  return 0;
}
