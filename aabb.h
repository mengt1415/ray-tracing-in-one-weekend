#ifndef AABB_H
#define AABB_H
#include "rtweekend.h"

class aabb
{
public:
	vec3 _min, _max;
	aabb(){}
	aabb(const point3& a, const point3& b)
	{
		_min= a;
		_max = b;
	}
	point3 min()const
	{
		return _min;
	}
	point3 max()const
	{
		return _max;
	}
	//t_exit>=t_enter && t_exit>0
	bool hit(const ray& r, double tmin, double tmax) const
	{
		for (int a = 0; a < 3; a++)
		{
			double invD = 1.0 / r.direction()[a];
			double t0 = (min()[a] - r.origin()[a]) * invD;
			double t1 = (max()[a] - r.origin()[a]) * invD;
			if (invD < 0.0)
				std::swap(t0, t1);
			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;
			if (tmax <= tmin)
				return false;
		}
		return true;
	}
	double area()
	{
		//double a = max().x() - min().x();
		double a = _max.x() - _min.x();
		double b = _max.y() - _min.y();
		double c = _max.z() - _min.z();
		return 2 * (a * a + b * b + c * c);
	}
	
	int longest_axis() const //in this const function,if vec3:: x() don't have const then error
	{
		double a = _max.x() - _min.x();
		double b = _max.y() - _min.y();
		double c = _max.z() - _min.z();
		if (a > b && a > c)
			return 0;
		else if (b > c)
			return 1;
		else
			return 2;

	}
};

aabb surrounding_box(aabb box0, aabb box1)
{
	vec3 small(fmin(box0.min().x(), box1.min().x()),
		fmin(box0.min().y(), box1.min().y()),
		fmin(box0.min().z(), box1.min().z()));

	vec3 big(fmax(box0.max().x(), box1.max().x()),
		fmax(box0.max().y(), box1.max().y()),
		fmax(box0.max().z(), box1.max().z()));

	return aabb(small, big);
}
#endif // !AABB_H

