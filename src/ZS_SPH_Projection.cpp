#include "ZS_Projection.h"

int main(int ,char** )
{
  time_t time_start,time_end;
  ZS_Projection sph_projection;
  time(&time_start);
  sph_projection.projection_init();
  sph_projection.projection_calculate();
  time(&time_end);
  cout << "Caculation Time : " << difftime(time_end,time_start) <<  endl;
  return EXIT_SUCCESS;
}
