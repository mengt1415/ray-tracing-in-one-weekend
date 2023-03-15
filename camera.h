#ifndef CAMERA_H
#define CAMERA_H
#include "ray.h"
double MY_PI = 3.1415926535897932;
class camera
{
public:
	camera(vec3 lookfrom,vec3 lookat,vec3 vup,double vfov,double aspect)
	{
		lens_radius = 0.0;
		dist_focus = 1.0;
		time0 = time1 = 0.0;
		double theta = vfov * MY_PI / 180;
		double half_height = tan(theta / 2);
		double half_width = aspect * half_height;

		origin = lookfrom;
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w,u);
		lower_left_corner = origin - half_height*v - half_width*u - w;
		vertical = 2 * half_height*v;
		horizental = 2 * half_width*u;
	}
	//with aperture dist_to_focus
	camera(vec3 lookfrom, vec3 lookat, vec3 vup, double vfov, double aspect, double aperture, double d_f)
	{
		lens_radius = aperture/2;
		dist_focus = d_f;
		double theta = vfov * MY_PI / 180;
		double half_height = tan(theta / 2);
		double half_width = aspect * half_height;

		origin = lookfrom;
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);
		lower_left_corner = origin - half_height *dist_focus* v - half_width*dist_focus * u - dist_focus*w;
		vertical = 2 * half_height*dist_focus * v;
		horizental = 2 * half_width*dist_focus * u;
	}

	// shutter can open during time interval
	camera(vec3 lookfrom, vec3 lookat, vec3 vup, double vfov, double aspect, double aperture, double d_f,double t0,double t1)
	{
		lens_radius = aperture / 2;
		dist_focus = d_f;
		time1 = t1;
		time0 = t0;
		double theta = vfov * MY_PI / 180;
		double half_height = tan(theta / 2);
		double half_width = aspect * half_height;

		origin = lookfrom;
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);
		lower_left_corner = origin - half_height * dist_focus * v - half_width * dist_focus * u - dist_focus * w;
		vertical = 2 * half_height * dist_focus * v;
		horizental = 2 * half_width * dist_focus * u;
	}

	ray get_ray(double s, double t)
	{
		vec3 rd = lens_radius * random_in_unit_disk();
		vec3 offset = u * rd.x() + v * rd.y();
		double time = time0 + my_drand48() * (time1 - time0);
		return ray(origin+offset, lower_left_corner + s* horizental + t * vertical - origin-offset,time);
	}
	vec3 u, v, w;
	point3 origin;
	point3 lower_left_corner;
	vec3 horizental;
	vec3 vertical;
	double lens_radius;
	double dist_focus;
	double time0;
	double time1;
};
#endif // !CAMERA_H
