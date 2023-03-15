#ifndef SPHERE_H
#define SPHERE_H

#include "rtweekend.h"
#include "hittable.h"

class sphere :public hittable
{
public:
	point3 cent;
	double radi;
	material* mat_ptr;
	sphere() {}
	sphere(point3 c, double r, material* m)
		:cent(c), radi(r), mat_ptr(m)
	{};
	virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
	virtual bool bounding_box(double t0, double t1, aabb& output_box) const override;
	vec3 center() const
	{
		return cent;
	}

	double radius()
	{
		return radi;
	}
	static void get_sphere_uv(const point3& p, double& u, double& v)
	{
		// p: a given point on the sphere of radius one, centered at the origin.
			// u: returned value [0,1] of angle around the Y axis from X=-1.
			// v: returned value [0,1] of angle from Y=-1 to Y=+1.
			//     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
			//     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
			//     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

		/*auto theta = acos(-p.y());
		auto phi = atan2(-p.z(), p.x()) + pi;

		u = phi / (2 * pi);
		v = theta / pi;*/
		double phi = atan2(p.z(), p.x());
		double theta = asin(p.y());
		u = 1 - (phi + pi) / (2 * pi);
		v = (theta + pi / 2) / pi;

	}
		

};

bool sphere::hit(const ray& r, double tmin, double tmax, hit_record& rec) const
{
	vec3 oc = r.origin() - cent;
	double a = dot(r.direction(), r.direction());
	double b = 2.0 * dot(oc, r.direction());
	double c = dot(oc, oc) - (radi) * (radi);
	float discriminant = b * b - 4* a * c;
	if (discriminant > 0)
	{
		// the first point the ray hit,that is the smaller root
		double temp = (-b - sqrt(discriminant)) / (2 * a);
		if (temp<tmax&&temp>tmin)
		{
			rec.t = temp;
			rec.p = r.at_t(rec.t);
			vec3 outward_normal = (rec.p - cent) / radi;
			rec.mat_ptr = mat_ptr;
			rec.normal = outward_normal;
			//rec.set_face_normal(r, outward_normal);
			get_sphere_uv(outward_normal, rec.u, rec.v);
			return true;
		}
		temp= (-b + sqrt(discriminant)) / (2 * a);
		if (temp<tmax && temp>tmin)
		{
			rec.t = temp;
			rec.p = r.at_t(rec.t);
			vec3 outward_normal = (rec.p - cent) / radi;
			rec.mat_ptr = mat_ptr;
			rec.normal = outward_normal;
			get_sphere_uv(outward_normal, rec.u, rec.v);
			//rec.set_face_normal(r, outward_normal);
			return true;
		}
	}
	return false;//without this will cause false image
	
}

bool sphere::bounding_box(double t0, double t1, aabb& output_box) const
{
	output_box = aabb(center() - vec3(radi, radi, radi), center() + vec3(radi, radi, radi));
	return true;
}

class m_sphere :public hittable
{
public:
	m_sphere(){}
	m_sphere(vec3 cen0, vec3 cen1, double r, material* m, double t0, double t1)
	{
		center0 = cen0;
		center1 = cen1;
		time0 = t0;
		time1 = t1;
		_radius = r;
		mat_ptr = m;
	}

	virtual bool bounding_box(double _time0, double _time1, aabb& output_box) const override{
		aabb box0(
			center(_time0) - vec3(_radius, _radius, _radius),
			center(_time0) + vec3(_radius, _radius, _radius));
		aabb box1(
			center(_time1) - vec3(_radius, _radius, _radius),
			center(_time1) + vec3(_radius, _radius, _radius));
		output_box = surrounding_box(box0, box1);
		return true;
	}

	virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const override //there must be const override here
	{																						//and const in center() otherwise error
		vec3 oc = r.origin() - center(r.time());
		double a = dot(r.direction(), r.direction());
		double b = dot(oc, r.direction());
		double c = dot(oc, oc) - _radius * _radius;
		double discri = b * b - a * c;
		if (discri > 0)
		{
			double temp = (-b - sqrt(discri)) / a;
			if (temp<tmax && temp>tmin)
			{
				rec.t = temp;
				rec.p = r.at_t(rec.t);
				rec.normal = (rec.p - center(r.time())) / _radius;
				rec.mat_ptr = mat_ptr;
				return true;
			}
			temp = (-b + sqrt(discri)) / a;
			if (temp<tmax && temp>tmin)
			{
				rec.t = temp;
				rec.p = r.at_t(rec.t);
				rec.normal = (rec.p - center(r.time())) / _radius;
				rec.mat_ptr = mat_ptr;
				return true;
			}
		}
		return false;
	}

	vec3 center(double time) const
	{
		return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
	}
	double radius()
	{
		return _radius;
	}
	vec3 center0;
	vec3 center1;
	double time0;
	double time1;
	double _radius;
	material* mat_ptr;

};

#endif // !SPHERE_H

