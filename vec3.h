#ifndef VEC3_H
#define VEC3_H
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<stdlib.h>
#include<iostream>
#include "rtweekend.h"
using std::sqrt;
using std::fabs;

//auto inline static
#define m_my_drand 0x100000000LL
#define c_my_drand 0xb16 // can't be 'c',for there are parameters in main.cpp called 'c'
#define a_my_drand 0x5deece66dll
static unsigned long long seed_my_drand = 1;
double my_drand48()
{
	
	//wrong!generate the same number,so MSAA fails
	//for (int i = 0; i < 10; i++)
		//cout << my_drand48() << endl;
	/*unsigned seed;
	seed = time(0);
	srand(seed);
	return rand() / double(RAND_MAX);*/
	seed_my_drand = (a_my_drand * seed_my_drand + c_my_drand) & 0xffffffffffffll;
	unsigned int x = seed_my_drand >> 16;
	return (double(x) / double(m_my_drand));


}

/*double random_double(double min, double max)
{
	unsigned seed;
	seed = time(0);
	srand(seed);
	return (rand() / double(RAND_MAX)) * (max - min) - min;
}*/

class vec3 {

public:
	double e[3];
	vec3(){}
	vec3(double e0, double e1, double e2)
	{
		e[0] = e0;
		e[1] = e1;
		e[2] = e2;
	}
	inline double x() const { return e[0]; }
	inline double y() const { return e[1]; }
	inline double z() const { return e[2]; }
	inline double r() const  { return e[0]; }
	inline double g() const  { return e[1]; }
	inline double b() const  { return e[2]; }
	
	vec3 operator-() const
	{
		return vec3(-e[0], -e[1], -e[2]);
	}

	double operator[](int i) const
	{
		return e[i];
	}

	double& operator[](int i)
	{
		return e[i];
	}
	vec3& operator+=(const vec3& v)
	{
		e[0] += v.e[0];
		e[1] += v.e[1];
		e[2] += v.e[2];
		return *this;
	}

	vec3& operator*=(const vec3& v)
	{
		e[0] *= v.e[0];
		e[1] *= v.e[1];
		e[2] *= v.e[2];
		return *this;
	}

	vec3& operator*=(const double t)
	{
		e[0] *= t;
		e[1] *= t;
		e[2] *= t;
		return *this;
	}

	vec3& operator/=(const double t)
	{
		e[0] /= t;
		e[1] /= t;
		e[2] /= t;
		return *this;
	}

	double length()
	{
		double s = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
		return sqrt(s);
	}

	double length_square()
	{
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}

	bool near_zero() const
	{
		const double s = 1e-8;
		return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
	}
	
	

	inline static vec3 random() {
		return vec3(random_double(), random_double(), random_double());
	}

	inline static vec3 random(double min, double max) {
		return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
	}
};

//#endif // !VEC3_H 会出问题 以下的函数C2048已有主体
using point3 = vec3;
using color = vec3;

inline std::ostream& operator<<(std::ostream& out, const vec3& v)
{
	return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3& u, const vec3& v) 
{
	return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3& u, const vec3& v) 
{
	return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3& u, const vec3& v) 
{
	return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3& v) 
{
	return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline vec3 operator*(const vec3& v, double t) 
{
	return t * v;
}

inline vec3 operator/(vec3 v, double t) 
{
	return (1 / t) * v;
}

inline double dot(const vec3& u, const vec3& v) 
{
	return u.e[0] * v.e[0]
		+ u.e[1] * v.e[1]
		+ u.e[2] * v.e[2];
}

inline vec3 cross(const vec3& u, const vec3& v)
{
	return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
		u.e[2] * v.e[0] - u.e[0] * v.e[2],
		u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
	return v / v.length();
}

inline vec3 random_in_unit_disk() {
	while (true) {
		auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
		if (p.length_square() >= 1) continue;
		return p;
	}
}

inline vec3 random_in_unit_sphere() {
	while (true) {
		auto p = vec3::random(-1, 1);
		if (p.length_square() >= 1) continue;
		return p;
	}
}

inline vec3 random_unit_vector() {
	return unit_vector(random_in_unit_sphere());
}

inline vec3 random_in_hemisphere(const vec3& normal) {
	vec3 in_unit_sphere = random_in_unit_sphere();
	if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
		return in_unit_sphere;
	else
		return -in_unit_sphere;
}

inline vec3 reflect(const vec3& v, const vec3& n) {
	return v - 2 * dot(v, n) * n;
}

inline vec3 refract1(const vec3& uv, const vec3& n, double etai_over_etat) {
	auto cos_theta = fmin(dot(-uv, n), 1.0);
	vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
	vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_square())) * n;//length_square
	return r_out_perp + r_out_parallel;
}

inline bool refract(const vec3& v, const vec3& n, double ni_over_nt, vec3& refreacted)
{
	vec3 uv = unit_vector(v);
	double dt = dot(uv, n);
	double discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
	if (discriminant > 0.0)
	{
		refreacted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
		return true;
	}
	else
		return false;
}

#endif // !VEC3_H