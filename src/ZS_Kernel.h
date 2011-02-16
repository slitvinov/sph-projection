#ifndef ZS_KERNEL_H
#define ZS_KERNEL_H
#include "ZS_Particles.h"

class ZS_Kernel : public ZS_Particles
{
	public:

		ZS_Kernel();
		~ZS_Kernel();

		double get_kernel(double r);
		double get_kernel_derivity(double r);
};

#endif //ZS_KERNEL_H
