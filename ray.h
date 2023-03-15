#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray
{
public:
	ray() {}
	ray(const point3& origin,const vec3& direction)
		:ori(origin),dir(unit_vector(direction)),_t(0)
	{}
	ray(const point3& origin,const vec3& direction,double time)
		:ori(origin),dir(direction),_t(time)
	{}
	point3 origin() const
	{
		return ori;
	}
	vec3 direction() const
	{
		return dir;
	}
	double time() const
	{
		return _t;
	}
	point3 at_t(double t) const
	{
		return ori + t * dir;
	}

	point3 ori;
	vec3 dir;
	double _t;
};

#endif
