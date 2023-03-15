#ifndef MATERIAL_H
#define MATERIAL_H
#include "rtweekend.h"
#include "hittable.h"
#include "texture.h"
struct hit_record;

class material
{
public:
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const = 0;
	virtual vec3 emitted(double u, double v, const vec3& p) const//
	{ //没有const导致渲染出来的图片全黑
		return vec3(0, 0, 0);
	}
};

class lambertian :public material
{
public:
	lambertian(texture *a) :albedo(a){}
	texture* albedo;
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered)const override
	{
		vec3 scatter_dir = rec.normal + random_unit_vector();
		//if (scatter_dir.near_zero())
			//scatter_dir = rec.normal;

		scattered = ray(rec.p, scatter_dir,r_in.time());
		attenuation = albedo->value(rec.u,rec.v,rec.p);
		return true;
	}
};
class metal :public material
{
public:
	color albedo;
	double fuzz;
	metal(const color& a,double f):albedo(a),fuzz(f<1?f:1){}
	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered)const override
	{
		vec3 reflect_dir = reflect(unit_vector(r_in.direction()), rec.normal);
		scattered = ray(rec.p, reflect_dir + fuzz * random_in_unit_sphere());
		attenuation = albedo;
		return(dot(scattered.direction(), rec.normal) > 0);
	}
};
//other DIELECTRIC
/*bool refract3(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted) {
	vec3 unit_v = unit_vector(v);
	float dt = dot(unit_v, n);//dt, or we say, cos theta
	float square_cos_theta1 = 1.0f - ni_over_nt * ni_over_nt * (1 - dt * dt);
	if (square_cos_theta1 > 0) {
		refracted = ni_over_nt * (unit_v - n * dt) - n * sqrt(square_cos_theta1);
		return true;
	}
	else
		return false;//
}

//Fresnel - Schlick approximation
float schlick(float cosine, float ref_idx) {
	float r0 = (1 - ref_idx) / (1 + ref_idx);
	r0 = r0 * r0;
	return r0 + (1 - r0) * pow((1 - cosine), 5);
}

class dielectric : public material {
public:
	dielectric(float ri) :ref_idx(ri) {

	}
	float ref_idx;
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered)const {
		vec3 outward_normal;
		vec3 reflected = reflect(r_in.direction(), rec.normal);
		float ni_over_nt;
		attenuation = vec3(1.0f, 1.0f, 1.0f);
		vec3 refracted;
		float reflect_prob;
		float cosine;

		if (dot(r_in.direction(), rec.normal) > 0) {
			outward_normal = -rec.normal;
			ni_over_nt = ref_idx;
			cosine = ref_idx * dot(r_in.direction(), rec.normal) / r_in.direction().length();
		}
		else {
			outward_normal = rec.normal;
			ni_over_nt = 1.0 / ref_idx;
			cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
		}

		if (refract3(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
			reflect_prob = schlick(cosine, ref_idx);
		}
		else {
			reflect_prob = 1.0f;
		}


		if (reflect_prob > (rand() / (float)(RAND_MAX))) {
			scattered = ray(rec.p, reflected);
		}
		else {
			scattered = ray(rec.p, refracted);
		}
		return true;
	}

};*/

//MY DIELECTRIC
 double schlick(double cos, double ref_idx)
{
	double r0 = (1 - ref_idx) / (1 + ref_idx);
	r0 = r0 * r0;
	return r0 + (1 - r0) * pow((1 - cos), 5);
}
class dielectric :public material
{
public:
	double ref_idx;
	dielectric(double r):ref_idx(r){}
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered)const override
	{
		vec3 outward_n;
		//vec3 reflected = reflect(r_in.direction(), rec.normal);
		vec3 reflect_dir = reflect(r_in.direction(), rec.normal);
		double ni_over_nt;
		double reflect_prob;
		attenuation = vec3(1.0, 1.0, 1.0);
		vec3 refract_dir;
		double cos;
		if (dot(r_in.direction(), rec.normal) > 0)
		{
			outward_n = -rec.normal;
			ni_over_nt = ref_idx;
			cos = ref_idx * dot(r_in.direction(), rec.normal) / r_in.direction().length();
		}
		else
		{
			outward_n = rec.normal;
			ni_over_nt = 1.0/ref_idx;
			cos = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
		}
		//reflect_dir = reflect(r_in.direction(), outward_n);
		if (refract(r_in.direction(), outward_n, ni_over_nt, refract_dir))//折射了，没有全反射
		{
			reflect_prob = schlick(cos,ref_idx);
		}
		else
		{
			reflect_prob = 1.0;
		}
		if (my_drand48() < reflect_prob)
		{
			scattered = ray(rec.p, reflect_dir);
		}
		else
		{
			scattered = ray(rec.p, refract_dir);
		}
		return true;

	}

	
};


class diffuse_light :public material
{
public:
	diffuse_light(texture* a):emit(a) {}
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const
	{
		return false;
	}
	virtual vec3 emitted(double u, double v, const vec3& p)const
	{
		return emit->value(u, v, p);
	}
	texture* emit;

};
//PeterShirly DIELECTRIC
/*class dielectric : public material {
public:
	dielectric(double index_of_refraction) : ir(index_of_refraction) {}

	virtual bool scatter(
		const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
	) const override {
		attenuation = color(1.0, 1.0, 1.0);
		double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

		vec3 unit_direction = unit_vector(r_in.direction());
		double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

		bool cannot_refract = refraction_ratio * sin_theta > 1.0;
		vec3 direction;

		if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
			direction = reflect(unit_direction, rec.normal);
		else
			direction = refract1(unit_direction, rec.normal, refraction_ratio);

		scattered = ray(rec.p, direction);
		return true;
	}

public:
	double ir; // Index of Refraction

private:
	static double reflectance(double cosine, double ref_idx) {
		// Use Schlick's approximation for reflectance.
		auto r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}
};*/

class isotropic :public material
{
public:
	isotropic(texture* a) :albedo(a) {}
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const
	{
		scattered = ray(rec.p, random_in_unit_sphere());
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}
	texture* albedo;
};

#endif // !MATERIAL_H
