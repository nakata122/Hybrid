// smallpaint by karoly zsolnai - zsolnai@cg.tuwien.ac.at
//
// render, modify, create new scenes, tinker around, and most of all:
// have fun!
//
// This program is used as an educational learning tool on the Rendering
// course at TU Wien. Course webpage:
// http://cg.tuwien.ac.at/courses/Rendering/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctime>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <cstdint>
#include <algorithm>
#include "../src_gui/main.h"

namespace smallpaint_fixed {

// Helpers for random number generation
std::mt19937 mersenneTwister;
std::uniform_real_distribution<double> uniform;

#define RND (2.0*uniform(mersenneTwister)-1.0)
#define RND2 (uniform(mersenneTwister))

#define PI 3.1415926536
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int width, height;
const double inf = 1e9;
const double eps = 1e-6;
using namespace std;
typedef unordered_map<string, double> pl;

struct Vec {
	double x, y, z;
	Vec(double x0 = 0, double y0 = 0, double z0 = 0) { x = x0; y = y0; z = z0; }
	Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x*b, y*b, z*b); }
	Vec operator/(double b) const { return Vec(x / b, y / b, z / b); }
	Vec mult(const Vec &b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec& norm() { return *this = *this * (1 / sqrt(x*x + y*y + z*z)); }
	double length() { return sqrt(x*x + y*y + z*z); }
	double dot(const Vec &b) const { return x*b.x + y*b.y + z*b.z; }
	Vec operator%(const Vec &b) const { return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }
	//	double& operator[](size_t i) { return data[i]; }
	const double& operator[](size_t i) const { return i == 0 ? x : (i == 1 ? y : z); }
};

// given v1, set v2 and v3 so they form an orthonormal system
// (we assume v1 is already normalized)
void ons(const Vec& v1, Vec& v2, Vec& v3) {
	if (std::abs(v1.x) > std::abs(v1.y)) {
		// project to the y = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.f / sqrtf(v1.x * v1.x + v1.z * v1.z);
		v2 = Vec(-v1.z * invLen, 0.0f, v1.x * invLen);
	} else {
		// project to the x = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.0f / sqrtf(v1.y * v1.y + v1.z * v1.z);
		v2 = Vec(0.0f, v1.z * invLen, -v1.y * invLen);
	}
	v3 = v1 % v2;
}

// Rays have origin and direction.
// The direction vector should always be normalized.
struct Ray {
	Vec o, d;
	Ray(Vec o0 = 0, Vec d0 = 0) { o = o0, d = d0.norm(); }
};

// Objects have color, emission, type (diffuse, specular, refractive)
// All object should be intersectable and should be able to compute their surface normals.
class Obj {
public:
	Vec cl;
	double emission;
	int type;
	void setMat(Vec cl_ = 0, double emission_ = 0, int type_ = 0) { cl = cl_; emission = emission_; type = type_; }
	virtual double intersect(const Ray&) const = 0;
	virtual Vec normal(const Vec&) const = 0;
};

class Plane : public Obj {
public:
	Vec n;
	double d;
	Plane(double d_ = 0, Vec n_ = 0) {
		d = d_;
		n = n_;
	}
	double intersect(const Ray& ray) const {
		double d0 = n.dot(ray.d);
		if (d0 != 0) {
			double t = -1 * (((n.dot(ray.o)) + d) / d0);
			return (t > eps) ? t : 0;
		} else return 0;
	}
	Vec normal(const Vec& p0) const { return n; }
};

class Sphere : public Obj {
public:
	Vec c;
	double r;

	Sphere(double r_ = 0, Vec c_ = 0) { c = c_; r = r_; }
    double intersect(const Ray& ray) const {
        Vec m = ray.o - c;
        float b = m.dot(ray.d);
        float c = m.dot(m) - r * r;

        // Exit if r’s origin outside s (c > 0) and r pointing away from s (b > 0)
        if (c > 0.0f && b > 0.0f) return 0;
        float discr = b*b - c;

        // A negative discriminant corresponds to ray missing sphere
        if (discr < 0.0f) return 0;

        // Ray now found to intersect sphere, compute smallest t value of intersection
        double t = -b - sqrt(discr);

        // If t is negative, ray started inside sphere so clamp t to zero
        if (t < 0.0f) t = 0.0f;

        return t;
	}

	Vec normal(const Vec& p0) const {
		return (p0 - c).norm();
	}
};

//returns a random point on the surface
Ray sphere(Obj* light, double u1, double u2) {
    Sphere* s = dynamic_cast<Sphere*>(light);
    double theta = 2 * PI * u1;
    double phi = u2 * PI;
    double x = sin(phi) * cos(theta) + s->c.x;
    double y = sin(phi) * sin(theta) + s->c.y;
    double z = cos(phi) + s->c.z;
    Vec d = Vec(x, y, z);
    return Ray(s->c, s->normal(d));
}

class Intersection {
public:
	Intersection() { t = inf; object = nullptr; }
	Intersection(double t_, Obj* object_) { t = t_; object = object_; }
	operator bool() { return object != nullptr; }
	double t;
	Obj* object;
};

class Scene {
	vector<Obj*> objects;

public:
	void add(Obj* object) {
		objects.push_back(object);
	}

	Intersection intersect(const Ray& ray) const {
		Intersection closestIntersection;
		// intersect all objects, one after the other
		for (auto iter = objects.begin(); iter != objects.end(); ++iter) {
			double t = (*iter)->intersect(ray);
			if (t > eps && t < closestIntersection.t) {
				closestIntersection.t = t;
				closestIntersection.object = *iter;
			}
		}
		return closestIntersection;
	}
};

// Class for generating the Halton low-discrepancy series for Quasi
// Monte Carlo integration.
class Halton {
	double value, inv_base;
public:
	void number(int i, int base) {
		double f = inv_base = 1.0 / base;
		value = 0.0;
		while (i > 0) {
			value += f * (double)(i%base);
			i /= base;
			f *= inv_base;
		}
	}
	void next() {
		double r = 1.0 - value - 0.0000001;
		if (inv_base < r) value += inv_base;
		else {
			double h = inv_base, hh;
			do { hh = h; h *= inv_base; } while (h >= r);
			value += hh + h - 1.0;
		}
	}
	double get() { return value; }
};

// Input is the pixel offset, output is the appropriate coordinate
// on the image plane
Vec camcr(const double x, const double y) {
	double w = width;
	double h = height;
	float fovx = PI / 4;
	float fovy = (h / w) * fovx;
	return Vec(((2 * x - w) / w) * tan(fovx),
			   -((2 * y - h) / h) * tan(fovy),
			   -1.0);
}

// Uniform sampling on a hemisphere to produce outgoing ray directions.
// courtesy of http://www.rorydriscoll.com/2009/01/07/better-sampling/
Vec hemisphere(double u1, double u2) {
	const double r = sqrt(1.0 - u1*u1);
	const double phi = 2 * PI * u2;
	return Vec(cos(phi)*r, sin(phi)*r, u1);
}

void trace(Ray &ray, const Scene& scene, int depth, Vec& flux, Vec **pix, pl& params) {
	// Russian roulette: starting at depth 5, each recursive step will stop with a probability of 0.1
	double rrFactor = 1.0;
    if (depth >= 5) {
        const double rrStopProbability = 0.1;
		if (RND2 <= rrStopProbability) {
			return;
		}
		rrFactor = 1.0 / (1.0 - rrStopProbability);
	}

	Intersection intersection = scene.intersect(ray);
	if (!intersection) return;


	// Travel the ray to the hit point where the closest object lies and compute the surface normal there.
    Vec hp = ray.o + ray.d * intersection.t;

    Vec N = intersection.object->normal(hp);
    ray.o = hp;
	// Add the emission, the L_e(x,w) part of the rendering equation, but scale it with the Russian Roulette
	// probability weight.
	const double emission = intersection.object->emission;
    flux = flux + Vec(emission, emission, emission) * rrFactor;


//    if (intersection.object->type == 4) {
//        Sphere* tp = dynamic_cast<Sphere*>(intersection.object);

//        if(hp.x >= 0 && hp.x < width && hp.y >= 0 && hp.y < height)
//            pix[(int)hp.x][(int)hp.y] = intersection.object->cl * 255;
//    }

	// Diffuse BRDF - choose an outgoing direction with hemisphere sampling.
    if (intersection.object->type == 1) {
        Vec rotX, rotY;
        ons(N, rotX, rotY);
        Vec sampledDir = hemisphere(RND2, RND2);
        Vec rotatedDir;
        rotatedDir.x = Vec(rotX.x, rotY.x, N.x).dot(sampledDir);
        rotatedDir.y = Vec(rotX.y, rotY.y, N.y).dot(sampledDir);
        rotatedDir.z = Vec(rotX.z, rotY.z, N.z).dot(sampledDir);
        ray.d = rotatedDir;	// already normalized
        double cost = ray.d.dot(N);
        Vec tmp;
        trace(ray, scene, depth + 1, flux, pix, params);
        flux = flux + (tmp.mult(intersection.object->cl)) * cost * 0.1 * rrFactor;
        float dist = 1 - (intersection.t);
        if(hp.x >= 0 && hp.x < width && hp.y >= 0 && hp.y < height)
            pix[(int)hp.x][(int)hp.y] = (intersection.object->cl.mult(flux)*(1. / PI)).mult(Vec(dist, dist, dist));
    }

//	// Specular BRDF - this is a singularity in the rendering equation that follows
//	// delta distribution, therefore we handle this case explicitly - one incoming
//	// direction -> one outgoing direction, that is, the perfect reflection direction.
//	if (intersection.object->type == 2) {
//		double cost = ray.d.dot(N);
//		ray.d = (ray.d - N*(cost * 2)).norm();
//		Vec tmp = Vec(0, 0, 0);
//		trace(ray, scene, depth + 1, tmp, params);
//        flux = flux + tmp * rrFactor;
//	}

//	// Glass/refractive BRDF - we use the vector version of Snell's law and Fresnel's law
//	// to compute the outgoing reflection and refraction directions and probability weights.
//	if (intersection.object->type == 3) {
//		double n = params["refr_index"];
//		double R0 = (1.0 - n) / (1.0 + n);
//		R0 = R0*R0;
//		if (N.dot(ray.d) > 0) { // we're inside the medium
//			N = N*-1;
//			n = 1 / n;
//		}
//		n = 1 / n;
//		double cost1 = (N.dot(ray.d))*-1; // cosine of theta_1
//		double cost2 = 1.0 - n*n*(1.0 - cost1*cost1); // cosine of theta_2
//		double Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0); // Schlick-approximation
//		if (cost2 > 0 && RND2 > Rprob) { // refraction direction
//			ray.d = ((ray.d*n) + (N*(n*cost1 - sqrt(cost2)))).norm();
//		} else { // reflection direction
//			ray.d = (ray.d + N*(cost1 * 2)).norm();
//		}
//		Vec tmp;
//		trace(ray, scene, depth + 1, tmp, params);
//        flux = flux + tmp * 1.15 * rrFactor;
//	}
}

void render(int id, int size, int spp, double refr_index) {

	srand(time(NULL));
	pl params;

	Scene scene;
	auto add = [&scene](Obj* s, Vec cl, double emission, int type) {
		s->setMat(cl, emission, type);
		scene.add(s);
	};

    // Radius, position, color, emission, type (1=diff, 2=spec, 3=refr) for spheres
    add(new Sphere(50, Vec(100, 100, -10)), Vec(1, 0, 0), 0, 1); // Middle sphere
    add(new Sphere(100, Vec(2.0, -2.05, -3.7)), Vec(1, 1, 1), 0, 1); // Right sphere
    //add(new Sphere(0.6, Vec(-1.75, -1.95, -3.1)), Vec(.01, .01, .5), 0, 1); // Left sphere
    // Position, normal, color, emission, type for planes
//    add(new Plane(2.5, Vec(0, 1, 0)), Vec(.1, .1, .1), 0, 1); // Bottom plane
//    add(new Plane(5.5, Vec(0, 0, 1)), Vec(.1, .1, .1), 0, 1); // Back plane
//    add(new Plane(2.75, Vec(1, 0, 0)), Vec(.5, .01, .01), 0, 1); // Left plane
//    add(new Plane(2.75, Vec(-1, 0, 0)), Vec(.01, .5, .01), 0, 1); // Right plane
//    add(new Plane(3.0, Vec(0, -1, 0)), Vec(.1, .1, .1), 0, 1); // Ceiling plane
    //add(new Plane(0, Vec(0, 0, -1)), Vec(1., 1., 0.), 0, 1); // Front plane

	params["refr_index"] = refr_index;
	params["spp"] = spp; // samples per pixel

	width = size;
	height = size;

	Vec **pix = new Vec*[width];
	for (int i = 0; i < width; i++) {
		pix[i] = new Vec[height];
	}

	// correlated Halton-sequence dimensions
	Halton hal, hal2;
	hal.number(0, 2);
	hal2.number(0, 2);

	bool running = true;

    Obj* light = new Sphere(1000, Vec(width/2, 0, -3));
    light->setMat(Vec(1, 1, 1), 50, 1);
    //scene.add(light);

    Ray ray;
    Vec flux = light->cl*light->emission*(PI*4.0);

    for (int s = 0; s < spp; s++) {
    //#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < 1000; i++) {
                if (running) {
                    ray = sphere(light, RND2, RND2);
                    Vec rotX, rotY;
                    ons(ray.d, rotX, rotY);
                    Vec sampledDir = hemisphere(RND2, RND2);
                    Vec rotatedDir;
                    rotatedDir.x = Vec(rotX.x, rotY.x, ray.d.x).dot(sampledDir);
                    rotatedDir.y = Vec(rotX.y, rotY.y, ray.d.y).dot(sampledDir);
                    rotatedDir.z = Vec(rotX.z, rotY.z, ray.d.z).dot(sampledDir);
                    ray.d = rotatedDir;
                    trace(ray, scene, 0, flux, pix, params);
                    if (!smallpaint::isRunning(id)) running = false;
                }
                /*if (running) {
                    Vec color;
                    Ray ray;
                    ray.o = (Vec(0, 0, 0)); // rays start out from here
                    Vec cam = camcr(col, row); // construct image plane coordinates
                    cam.x = cam.x + RND / 700; // anti-aliasing for free
                    cam.y = cam.y + RND / 700;
                    ray.d = (cam - ray.o).norm(); // point from the origin to the camera plane
                    trace(ray, scene, 0, color, params);
                    pix[col][row] = pix[col][row] + color;
                    if (!smallpaint::isRunning(id)) running = false;
                }*/
        }
        if (!running) return;
        imageOutput(pix, s);
    }

}
}
