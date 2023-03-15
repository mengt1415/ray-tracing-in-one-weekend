#ifndef TRANSLATE_H
#define TRANSLATE_H
#include "rtweekend.h"
#include "hittable.h"
//translate and rotate
class translate :public hittable
{
public:
	translate(hittable* p, const vec3& displacement)
		:ptr(p), offset(displacement)
	{}
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;
	virtual bool bounding_box(double t0, double t1, aabb& box) const;
	hittable* ptr;
	vec3 offset;
};

bool translate::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
	ray moved_r(r.origin() - offset, r.direction(), r.time());
	if (ptr->hit(moved_r, t_min, t_max, rec))
	{
		rec.p += offset;
		return true;
	}
	else
		return false;
}

bool translate::bounding_box(double t0, double t1, aabb& box) const
{
	if (ptr->bounding_box(t0, t1, box))
	{
		box = (aabb(box.min() + offset, box.max() + offset));
		return true;
	}
	else
		return false;
}

class rotate_y :public hittable
{
public:
	rotate_y(hittable* p, double angle);
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;
	virtual bool bounding_box(double t0, double t1, aabb& _box) const
	{
		_box = box;
		return hasbox;
	}
	hittable* ptr;
	double sin_t;
	double cos_t;
	bool hasbox;
	aabb box;
};

rotate_y::rotate_y(hittable* p, double angle)
	:ptr(p)
{
	double radians = (pi / 180.) * angle;
	sin_t = sin(radians);
	cos_t = cos(radians);
	hasbox = ptr->bounding_box(0, 1, box);
	vec3 min(DBL_MAX, DBL_MAX, DBL_MAX);
	vec3 max(-DBL_MAX, -DBL_MAX, -DBL_MAX);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				double x = i * box.max().x() + (1 - i) * box.min().x();
				double y = j * box.max().y() + (1 - j) * box.min().y();
				double z = k * box.max().z() + (1 - k) * box.min().z();
				double newx = cos_t * x + sin_t * z;
				double newz = -sin_t * x + cos_t * z;
				vec3 tester(newx, y, newz);
				for (int c = 0; c < 3; c++)
				{
					if (tester[c] > max[c])
						max[c] = tester[c];
					if (tester[c] < min[c])
						min[c] = tester[c];
				}
			}
		}
	}
	box = aabb(min, max);
}

bool rotate_y::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
	vec3 origin = r.origin();
	vec3 direction = r.direction();
	origin[0] = cos_t * r.origin()[0] - sin_t * r.origin()[2];
	origin[2] = sin_t * r.origin()[0] + cos_t * r.origin()[2];
	direction[0] = cos_t * r.direction()[0] - sin_t * r.direction()[2];
	direction[2] = sin_t * r.direction()[0] + cos_t * r.direction()[2];
	ray rotated_r(origin, direction, r.time());
	if (ptr->hit(rotated_r, t_min, t_max, rec))
	{
		vec3 p = rec.p;
		vec3 normal = rec.normal;
		p[0] = cos_t * rec.p[0] + sin_t * rec.p[2];
		p[2] = -sin_t * rec.p[0] + cos_t * rec.p[2];
		normal[0] = cos_t * rec.normal[0] + sin_t * rec.normal[2];
		normal[2] = -sin_t * rec.normal[0] + cos_t * rec.normal[2];
		rec.p = p;
		rec.normal = normal;
		return true;
	}
	else
		return false;

}


#endif // !TRANSLATE_H

