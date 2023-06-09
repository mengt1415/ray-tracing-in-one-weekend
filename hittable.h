#ifndef HITTABLE_H
#define HITTABLE_H

#include "rtweekend.h"
//#include "hittable_list.h"
//本来将rect box都写在这里，这样就写box就需要include hittable——list，可是
//这样就报错了
//hittablelist那边有include hittable。h
#include "aabb.h"

class material;

struct hit_record
{
	point3 p;
	vec3 normal;
	material* mat_ptr;
	double u;
	double v;
	double t;

	bool front_face;

	inline void set_face_normal(const ray& r, const vec3& outward_normal)
	{
		front_face = dot(r.direction(), outward_normal) < 0;
		normal = front_face ? outward_normal : -outward_normal;
	}
};

class hittable
{
public:
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
	virtual bool bounding_box(double t0, double t1, aabb& box) const = 0;
};


#endif