#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include "time.h"
//#include<stdlib.h>
#include<thread>
#include<stdio.h>
#include<cmath>
#include<vector>
#include "rtweekend.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include <omp.h>
#include "bvh.h"
#include "rectbox.h"
#include "constant_medium.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "translate.h"
//#include "vec3.h"
//#include "ray.h"
// if keep both "#include vec3 ray" in main.cpp and rtweekend.h then errors
using namespace std;
double MAXDOUBLE = 999999999.99;
/*class my_sphere
{
public:
	//sphere(){}
	my_sphere(vec3 c, double r,vec3 col)
		:center(c), radius(r),color(col)
	{}
	point3 center;
	double radius;
	vec3 color;
	//sphere() :center(vec3(1.0, 1.0, 1.0)), radius(1.0) {}
	
	
};
my_sphere* slist[11];
my_sphere a(vec3(0.0, 0.0, -1.0), 0.25, vec3(1.0, 0.0, 0.0)),
b(vec3(0.25, 0, -2.0),0.5, vec3(0.0, 1.0, 0.0)),
c(vec3(0.0,0.25, -3),1.5, vec3(0.0, 0.0, 1.0)),
d(vec3(-3,0,-3),1.0,vec3(1.0,0.0,1.0)),
e(vec3(-6, 0, -10), 1.0, vec3(1.0, 0.0, 1.0)),
f(vec3(-9, 0, -10), 1.0, vec3(1.0, 0.0, 1.0)),
g(vec3(-9, 0, -10), 1.0, vec3(1.0, 0.0, 1.0)),
h(vec3(-12, 0, -10), 1.0, vec3(1.0, 0.0, 1.0)),
i(vec3(-15, 0, -10), 1.0, vec3(1.0, 0.0, 1.0)),
j(vec3(-18, 0, -10), 1.0, vec3(1.0, 0.0, 1.0)),
k(vec3(-21, 0, -10), 1.0, vec3(1.0, 0.0, 1.0));*/

/*vec3 r_color(ray& r)
{
	vec3 unit_direction = unit_vector(r.direction());
	double t = 0.5 * (unit_direction.y() + 1);
	return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 0.1);
}*/
/*bool hit_sphere(my_sphere* s[], ray& r, my_sphere*& hit_s)
{
	double tmin = 1e9;
	my_sphere* sp = nullptr;
	for (int i = 0; i <=3; i++)
	{
		vec3 oc = r.origin() - s[i]->center;
		double a = dot(r.direction(), r.direction());
		double b = 2.0 * dot(oc, r.direction());
		double c = dot(oc, oc) - (s[i]->radius) * (s[i]->radius);
		float discriminant = b * b - 4 * a * c;
		if (discriminant > 0)
		{
			double t = (-b + sqrt(discriminant)) / (2 * a);
			if (t < tmin)
			{
				r.t = t;
				tmin = t;
				sp = s[i];
				hit_s = s[i];
			}
		}
	}
	if (sp == nullptr)
		return false;
	else
		return true;


}*/
/*vec3 my_r_color(ray& r)
{
	my_sphere* s;
	if (hit_sphere(slist,r,s))
	{
		//return s->color;
		vec3 n = unit_vector(r.at_t(r.t) - s->center);
		return 0.5 * vec3(n.x() + 1, n.y() + 1, n.z() + 1);
	}
	vec3 unit_dir = unit_vector(r.direction());
	double t = 0.5 * (unit_dir.y() + 1);
	return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
}*/

vec3 r_color(const ray& r, hittable* world,int depth)
{   //hittable->materal decides how rays are scattered
	hit_record rec;
	
	if (world->hit(r, 0.001, MAXDOUBLE, rec)) //0.001 is important book1 chap7 the picture got much more brighter
	{
		//模拟光线随机的散射，当然每次散射一根光线，每个像素有大的w
		ray scattered;
		vec3 attenuation;
		vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
		if (depth < 16 && rec.mat_ptr->scatter(r, rec, attenuation, scattered))
		{
			
			return emitted+attenuation * r_color(scattered, world, depth + 1);
		}
		else
			return emitted;
	}
	else{return vec3(0, 0, 0);}
	//注释掉else部分，有了黑暗的背景
	/*else
	{
		vec3 unit_dir = unit_vector(r.direction());
		double t = 0.5 * (unit_dir.y() + 1.0);
		return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
	}*/
}
void threadmain(int thread_num,int thread_id,vector<vec3>& framebuffer, hittable* world,camera &mycamera,const int &nx,const int &ny,const int &spp)
{
	for (int j = ny - thread_id-1; j >= 0; j-=thread_num)
	{
		for (int i = 0; i < nx; i++)
		{

			vec3 col(0.0, 0.0, 0.0);
			for (int k = 0; k < spp; k++)
			{   //drand48() linux
				double u = double(i + my_drand48()) / double(nx); //x m dx
				double v = double(j + my_drand48()) / double(ny); //y n dy
				ray r = mycamera.get_ray(u, v);
				col += r_color(r, world, 0);
			}
			col /= double(spp);
			col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
			framebuffer[(ny - 1 - j) * nx + i] = col;
			//float r = float(i) / float(nx);
				//float g = float(j) / float(ny);
				//float b = 0.2;
				//vec3 col=my_r_color(r);
				//vec3 col = r_color(r,world);
				//int ir = int(255.99 * col[0]);
				//int ig = int(255.99 * col[1]);
				//int ib = int(255.99 * col[2]);
				//oppm << ir << " " << ig << " " << ib << "\n";
		}
	}
}
void mtrt(vector<vec3>& framebuffer, hittable* world, camera& mycamera, const int& nx, const int& ny, const int& spp)
{
	const int threadnum = thread::hardware_concurrency() - 2;
	vector<thread> threads;
	cout << "rt using " << threadnum << " threads\n";
	for (int i = 0; i != threadnum; i++)
	{
		threads.emplace_back(thread{ threadmain,threadnum,i,ref(framebuffer),world,ref(mycamera),cref(nx),cref(ny),cref(spp) });
	}
	for (auto& thread : threads)
	{
		thread.join();
	}
}
void omprt(vector<vec3>&framebuffer,hittable *world,camera &mycamera,int nx,int ny,int spp)
{
	int core = omp_get_num_procs();
	cout<<"ompRT omp_get_num_procs:" << core << endl;
	#pragma omp parallel for

	for (int j = ny - 1; j >= 0; j--)   
	{
		for (int i = 0; i < nx; i++)
		{

			vec3 col(0.0, 0.0, 0.0);
			for (int k = 0; k < spp; k++)
			{   //drand48() linux
				double u = double(i + my_drand48()) / double(nx); //x m dx
				double v = double(j + my_drand48()) / double(ny); //y n dy
				ray r = mycamera.get_ray(u, v);
				col += r_color(r, world, 0);
			}
			col /= double(spp);
			col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
			framebuffer[(ny - 1 - j) * nx + i] = col;
			//float r = float(i) / float(nx);
				//float g = float(j) / float(ny);
				//float b = 0.2;
				//vec3 col=my_r_color(r);
				//vec3 col = r_color(r,world);
				//int ir = int(255.99 * col[0]);
				//int ig = int(255.99 * col[1]);
				//int ib = int(255.99 * col[2]);
				//oppm << ir << " " << ig << " " << ib << "\n";
		}
	}
}
void output_pixel(const vector<vec3>& framebuffer,const int&nx,const int &ny)
{
	FILE* fp = fopen("fpic-camera.ppm", "wb");
	(void)fprintf(fp, "P6\n%d %d\n255\n", nx, ny);
	#pragma omp parallel for
	for (int i = 0; i < nx * ny; i++)
	{
		static unsigned char color[3];
		//没有clamp 会出现离谱的结果！
		color[0] = (unsigned char)(255.99 * clamp(framebuffer[i].x(), 0.0, 0.999));
		color[1] = (unsigned char)(255.99 * clamp(framebuffer[i].y(), 0.0, 0.999));
		color[2] = (unsigned char)(255.99 * clamp(framebuffer[i].z(), 0.0, 0.999));
		fwrite(color, 1, 3, fp);
	}
	fclose(fp);
}

hittable_list* random_scene()
{
	int n = 500;
	int nx, ny, nn;
	unsigned char* tex_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);
	material* mat = new lambertian(new image_texture(tex_data, nx, ny));
	texture* pertext = new noise_texture(5);
	hittable** list = new hittable * [n + 1];
	texture* checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), new constant_texture(vec3(0.9, 0.9, 0.9)));
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(checker));
	int i = 1;
	for (int a = -11; a < 11; a++)
	{
		for (int b = -11; b < 11; b++)
		{



			double choose_mat = my_drand48();
			vec3 center(a + 0.9 * my_drand48(), 0.2, b + 0.9 * my_drand48());
			if ((center - vec3(4, 0.2, 0)).length()>0.9)
			{
				if (choose_mat < 0.1)
				{
					list[i++] = new sphere(center, 0.2, new lambertian(pertext));
				}
			else if(choose_mat<0.6)
				{
					//diffuse
					list[i++] = new sphere(center, 0.2, new lambertian(new constant_texture(vec3(my_drand48() * my_drand48(), my_drand48() * my_drand48(), my_drand48() * my_drand48()))));

				}
			else if (choose_mat < 0.75)
			{
				//metal
				list[i++] = new sphere(center, 0.2, new metal(vec3(0.5 * (1 + my_drand48()), 0.5 * (1 + my_drand48()), 0.5 * (1 + my_drand48())), 0.5 * my_drand48()));
			}
			else
			{
				//glass
				list[i++] = new sphere(center, 0.2, new dielectric(1.5));
			}
			}
		}
	}
	list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
	//list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
	list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(pertext));
	//list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
	list[i++] = new sphere(vec3(4, 1, 0), 1.0, mat);
	return new hittable_list(list, i);
}
//motion  blur
hittable_list* m_random_scene()
{
	int n = 500;
	hittable** list = new hittable * [n + 1];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(0.5, 0.5, 0.5))));
	int i = 1;
	for (int a = -10; a < 10; a++)
	{
		for (int b = -10; b < 10; b++)
		{
			double choose_mat = my_drand48();
			vec3 center(a + 0.9 * my_drand48(), 0.2, b + 0.9 * my_drand48());
			if ((center - vec3(4, 0.2, 0)).length() > 0.9)
			{
				if (choose_mat < 0.7)
				{
					//diffuse moving_sphere
					list[i++] = new m_sphere(center,center+vec3(0,0.5*my_drand48(),0), 0.2, new lambertian(new constant_texture(vec3(my_drand48() * my_drand48(), my_drand48() * my_drand48(), my_drand48() * my_drand48()))),1.0,2.0);

				}
				else if (choose_mat < 0.9)
				{
					//metal
					list[i++] = new sphere(center, 0.2, new metal(vec3(0.5 * (1 + my_drand48()), 0.5 * (1 + my_drand48()), 0.5 * (1 + my_drand48())), 0.5 * my_drand48()));
				}
				else
				{
					//glass
					list[i++] = new sphere(center, 0.2, new dielectric(1.5));
				}
			}
		}
	}
	list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
	list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(new constant_texture(vec3(0.4, 0.2, 0.1))));
	list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
	return new hittable_list(list, i);
}
//perlin noise
hittable_list* two_perlin_sphere()
{
	texture* pertext = new noise_texture(5);
	hittable** list = new hittable * [2];
	int nx, ny, nn;
	unsigned char* tex_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);
	material* mat = new lambertian(new image_texture(tex_data, nx, ny));
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
	return new hittable_list(list, 2);
}
hittable_list* image_texture_sphere()
{
	texture* pertext = new noise_texture(5);
	hittable** list = new hittable * [3];
	int nx, ny, nn;
	unsigned char* tex_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);
	material* mat = new lambertian(new image_texture(tex_data, nx, ny));
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new sphere(vec3(2, 1, 0), 1, mat);
	list[2] = new sphere(vec3(-2, 1, 0), 1,new metal(vec3(0.7, 0.6, 0.5), 0.0));
	return new hittable_list(list, 2);
}

hittable_list* simple_light()
{
	texture* pertext = new noise_texture(5);
	hittable** list = new hittable * [3];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
	list[2] = new xy_rect(3, 5, 1, 3, -2, new diffuse_light(new constant_texture(vec3(4, 4, 4))));
	return new hittable_list(list, 3);
}

hittable_list* cornell_box()
{
	hittable** list = new hittable* [20];
	int i = 0;
	material* red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
	material* white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
	material* green = new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
	material* light = new diffuse_light(new constant_texture(vec3(15, 15, 15)));
	//material* light = new diffuse_light(new constant_texture(vec3(15, 15, 15)));
	list[i++] = new flip_normal(new yz_rect(0, 555, 0, 555, 555, green));//0
	list[i++] = new yz_rect(0, 555, 0, 555, 0, red);//1
	list[i++] = new xz_rect(213, 343, 227, 332, 554, light);//2 +-1
	list[i++] = new flip_normal(new xz_rect(0, 555, 0, 555, 555, white));//3
	list[i++] = new xz_rect(0, 555, 0, 555, 0, white);//4
	list[i++] = new flip_normal(new xy_rect(0, 555, 0, 555, 555, white));//5
	//list[i++] = new box(vec3(130, 0, 65), vec3(295, 165, 230),white );//6
	//list[i++] = new box(vec3(265, 0, 295), vec3(430, 330, 460), white);//7
	list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 165, 165), white), -18), vec3(130, 0, 65));
	list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 330, 165), white), 15), vec3(265, 0, 295));
	return new hittable_list(list, i);
}
hittable_list* cornell_smoke()
{
	hittable** list = new hittable * [20];
	int i = 0;
	material* red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
	material* white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
	material* green = new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
	material* light = new diffuse_light(new constant_texture(vec3(15, 15, 15)));
	list[i++] = new flip_normal(new yz_rect(0, 555, 0, 555, 555, green));//0
	list[i++] = new yz_rect(0, 555, 0, 555, 0, red);//1
	list[i++] = new xz_rect(113, 443, 127, 432, 554, light);//2
	list[i++] = new flip_normal(new xz_rect(0, 555, 0, 555, 555, white));//3
	list[i++] = new xz_rect(0, 555, 0, 555, 0, white);//4
	list[i++] = new flip_normal(new xy_rect(0, 555, 0, 555, 555, white));//5
	//list[i++] = new box(vec3(130, 0, 65), vec3(295, 165, 230),white );//6
	//list[i++] = new box(vec3(265, 0, 295), vec3(430, 330, 460), white);//7
	list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 165, 165), new metal(vec3(0.8,0.8,0.9),0.1)), -18), vec3(130, 0, 65));
	//hittable *b2=
	hittable *b2 = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 330, 165), new metal(vec3(0.8, 0.8, 0.9), 0.01)), 15), vec3(265, 0, 295));
	list[i++] = b2;
	//list[i++] = new constant_medium(b2, 0.01, new constant_texture(vec3(0.0, 0.0, 0.0)));
	return new hittable_list(list, i);
}

hittable_list* final()
{
	int nb = 20;
	hittable** list = new hittable * [30];
	hittable** boxlist = new hittable * [10000];
	hittable** boxlist2 = new hittable * [10000];
	material* white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
	material* ground = new lambertian(new constant_texture(vec3(0.48, 0.83, 0.53)));
	int b = 0;
	for (int i = 0; i < nb; i++)
	{
		for (int j = 0; j < nb; j++)
		{
			double w = 100;
			double x0 = -1000 + i * w;
			double z0 = -1000 + j * w;
			double y0 = 0;
			double x1 = x0 + w;
			double y1 = 100 * (my_drand48() + 0.01);
			double z1 = z0 + w;
			boxlist[b++] = new box(vec3(x0, y0, z0), vec3(x1, y1, z1), ground);
		}
	}
	int l = 0;
	list[l++] = new bvh_node(boxlist, b, 0, 1);

	material* light = new diffuse_light(new constant_texture(vec3(7, 7, 7)));

	list[l++] = new xz_rect(123, 423, 147, 412, 554, light);

	vec3 center(400, 400, 200);

	list[l++] = new m_sphere(center, center + vec3(30, 0, 0),50, new lambertian(new constant_texture(vec3(0.57,0.44,0.86))),0,1);
	
	list[l++] = new sphere(vec3(260, 150, 45), 50, new dielectric(1.5));
	
	list[l++] = new sphere(vec3(0, 150, 145), 50, new metal(vec3(0.8, 0.8, 0.9), 0.02));
	
	hittable* boundary = new sphere(vec3(360, 150, 145), 70, new dielectric(1.5));
	list[l++] = boundary;
	
	list[l++] = new constant_medium(boundary, 0.2, new constant_texture(vec3(0.2, 0.4, 0.9)));
	
	boundary = new sphere(vec3(0, 0, 0), 5000, new dielectric(1.5));
	
	list[l++] = new constant_medium(boundary, 0.0001, new constant_texture(vec3(1.0, 1.0, 1.0)));
	
	int nx, ny, nn;
	
	unsigned char* tex_data = stbi_load("earthmap.jpg", &nx, &ny, &nn,0);
	
	material* emat = new lambertian(new image_texture(tex_data, nx, ny));
	
	list[l++] = new sphere(vec3(400, 200, 400), 100, emat);
	texture* pertext = new noise_texture(10);
	list[l++] = new sphere(vec3(220, 280, 300), 80, new lambertian(pertext));
	int ns = 1000;
	for (int j = 0; j < ns; j++)
	{
		boxlist2[j] = new sphere(vec3(165 * my_drand48(), 165 * my_drand48(), 165 * my_drand48()), 10, white);
		
	}
	list[l++] = new translate(new rotate_y(new bvh_node(boxlist2, ns, 0.0, 1.0), 15), vec3(-100, 270, 395));
	return new hittable_list(list, l);

}
int main()
{

	//ofstream oppm("pic-1.ppm", ios::out);
	//hittable* list[5];
	//int sx, sy, ss;
	//unsigned char* tex_data = stbi_load("earthmap.jpg", &sx, &sy, &ss, 0);
	//material* mat = new lambertian(new image_texture(tex_data, sx, sy));
	//list[0] = new sphere(vec3(0, 0, -1), 0.5,new lambertian(vec3(0.8,0.3,0.3)) ); wrong with shared_ptr? change to material* 
	//list[0] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(new constant_texture(vec3(0.8, 0.3, 0.3))));
	//list[0] = new sphere(vec3(0, 0, -1.5), 0.5, new metal(vec3(0.7,0.6,0.5),0.0));
	//list[1] = new sphere(vec3(0, -100.5, -1), 100.0, new lambertian(new constant_texture(vec3(0.8, 0.8, 0.0))));
	//list[2] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.1));
	//list[2] = new sphere(vec3(1, 0, -1), 0.5, mat);
	//list[3] = new sphere(vec3(-1, 0, -1), 0.5, new dielectric(1.5));
	//list[4] = new sphere(vec3(-1, 0, -1), -0.45, new dielectric(1.5));
	//hittable_list* world = new hittable_list(list, 5);
	clock_t start, end;
	start = clock();
	int nx = 1000, ny = 1000;
	//int nx = 2560,ny = 1600;
	
	int spp = 20000;
	//oppm << "P3\n" << nx << " " << ny << "\n255\n";
	point3 lower_left_cornetr(-2.0, -1.0, -1.0);
	vec3 horizontal(4.0, 0.0, 0.0);
	vec3 vertical(0.0, 2.0, 0.0);
	vec3 origin(0, 0, 0);

	/*vec3 lookfrom(0, 0, 2);
	vec3 lookat(0, 0, -1);
	vec3 vup(0, 1, 0);
	double vfov=20.0;
	camera mycamera(lookfrom, lookat, vup, 30, nx / (double)ny);*/

	/*vec3 lookfrom(3, 3, 2);
	vec3 lookat(0, 0, -1);
	vec3 vup(0, 1, 0);
	double vfov=20.0;
	double dist_focus = (lookfrom - lookat).length();
	double aperture = 0.2;
	camera mycamera(lookfrom, lookat, vup, 20,double(nx) / (double)ny, 0.0, 1);*/

	//random scene two perlin_niose sphere
	//vec3 lookfrom(13, 2, 3);
	//vec3 lookat(0, 0, 0);
	// double vfov=20.0;
	//vec3 vup(0, 1, 0);

	//simple_light
	/*vec3 lookfrom(26, 3, 6);
	nx = 5120, ny = 3200;
	vec3 lookat(0, 2, 0);
	vec3 vup(0, 1, 0);
	double vfov = 20.0;

	double dist_focus = 10.0;
	double aperture = 0.0;*/
	//camera mycamera(lookfrom, lookat, vup, vfov, double(nx) / (double)ny, aperture, dist_focus,1.0,2.0);
	
	//cornell box
	/*nx = 200, ny = 200;
	vec3 lookfrom(278, 278, -700);
	vec3 lookat(278, 278, 0);
	vec3 vup(0, 1, 0);
	double vfov = 40.0;
	double dist_focus = 10.0;
	double aperture = 0.02;*/

	//final_scene
	nx = 5120, ny = 3200;
	vec3 lookfrom(478, 278, -600);
	vec3 lookat(278, 278, 0);
	double vfov = 40.0;
	vec3 vup(0, 1, 0);
	double aperture = 0.0;
	double dist_focus = 10.0;
	camera mycamera(lookfrom, lookat, vup, vfov, double(nx) / (double)ny, aperture, dist_focus, 1.0, 2.0);
	hittable_list* world = final();
	vector<vec3> framebuffer(nx * ny);
	mtrt(framebuffer, world, mycamera, nx, ny, spp);
	//omprt(framebuffer, world, mycamera, nx, ny, spp);

	/*slist[0] = &a;
	slist[1] = &b;
	slist[2] = &c;
	slist[3] = &d;
	slist[4] = &e;
	slist[5] = &f;
	slist[6] = &g;
	slist[7] = &h;
	slist[8] = &i;
	slist[9] = &j;
	slist[10] = &k;*/
	/*int core = omp_get_num_procs();
	cout << core << endl;
    #pragma omp parallel for
	
	for (int j = ny - 1; j >= 0; j--)
	{
		for (int i = 0; i < nx; i++)
		{

			vec3 col(0.0, 0.0, 0.0);  
			for (int k = 0; k < spp; k++)
			{   //drand48() linux
				double u = double(i + my_drand48()) / double(nx); //x m dx
				double v = double(j + my_drand48()) / double(ny); //y n dy
				ray r = mycamera.get_ray(u, v);
				col += r_color(r, world, 0);
			}
			col /= double(spp);
			col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
			framebuffer[(ny-1-j) * nx + i] = col;
			//float r = float(i) / float(nx);
				//float g = float(j) / float(ny);
				//float b = 0.2;
				//vec3 col=my_r_color(r);
				//vec3 col = r_color(r,world);
				//int ir = int(255.99 * col[0]);
				//int ig = int(255.99 * col[1]);
				//int ib = int(255.99 * col[2]);
				//oppm << ir << " " << ig << " " << ib << "\n";
		}
	}*/
	output_pixel(framebuffer, nx, ny);
	/*FILE* fp = fopen("fpic-camera.ppm", "wb");
	(void)fprintf(fp, "P6\n%d %d\n255\n", nx, ny);
	#pragma omp parallel for
	for (int i=0;i<nx*ny; i++)
	{
		static unsigned char color[3];
		//没有clamp 会出现离谱的结果！
		color[0] = (unsigned char)(255.99 * clamp(framebuffer[i].x(),0.0,0.999));
		color[1] = (unsigned char)(255.99 * clamp(framebuffer[i].y(),0.0,0.999));
		color[2] = (unsigned char)(255.99 * clamp(framebuffer[i].z(),0.0,0.999));
		fwrite(color, 1, 3, fp);
	}
	fclose(fp);*/
		//for (int i = 0; i <= 100; i++)
		//{
		//	cout<<my_drand48()<<endl;
		//}
	/*my_sphere* s[2];
	my_sphere a, b(vec3(2.0, 0, 0), 2.0);
	s[0] = &a, s[1] = &b;
	cout << s[0]->center << " " << s[1]->center<<endl;// my_sphere* s[2], can not s[0].center
	s[0]->center = vec3(5, 5, 5);
	cout << s[0]->center << " " << a.center;*/
	//cout << RAND_MAX;
	/*my_sphere* cc = &a;
	cout << cc->center;
	test(cc);
	cout << endl << cc->center;*/
	//for (int i = 0; i < 10000; i++)
		//cout << my_drand48() << endl;
	end = clock();
	cout << "\n time: " << end - start << "\n";
	return 0;
}