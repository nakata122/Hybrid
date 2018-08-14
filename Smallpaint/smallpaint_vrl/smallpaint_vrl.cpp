// smallpaint by karoly zsolnai - zsolnai@cg.tuwien.ac.at
// virtual ray lights implementation by michael oppitz - e1227129@student.tuwien.ac.at
//
// g++ smallpaint_vrl.cpp -O3 -std=gnu++0x -fopenmp -static-libgcc -static-libstdc++
//
// This is an implementation of the "Virtual Ray Lights for Rendering Scenes with Participating Media"
// paper by Nov�k et al.
//
// This program is used as an educational learning tool on the Rendering
// course at TU Wien. Course webpage:
// http://cg.tuwien.ac.at/courses/Rendering/

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <assert.h>
#include <limits>
#include <iostream>
#include <fstream>
#include "../src_gui/main.h"

namespace smallpaint_vrl {

// Helpers for random number generation
std::mt19937_64 mersenneTwister;
std::uniform_real_distribution<double> uniform;

#define RND (2.0*uniform(mersenneTwister)-1.0)
#define RND2 (uniform(mersenneTwister))

#define PI 3.1415926536
#define INV_4PI 1.0/(4.0*PI)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define DIFFUSE 1
#define SPECULAR 2
#define REFRACTIVE 3

#define LIGHTPASS 0
#define RENDERPASS 1

#define LIGHT_BIDIR 0
#define LIGHT_DIR 1

#define SAMPLE_UNI 0
#define SAMPLE_EXP 1
#define SAMPLE_ORJ 2
#define SAMPLE_SIJ 3
#define SAMPLE_ADJ 4

int width, height;
const double inf = 1e3;
const double eps = 1e-6;
using namespace std;

struct RenderSettings {
	int scene = -1;
	int spp = 100, vrlps = 1;
	double sigma_a = -1, sigma_s = -1, g = -1;
	int bounces = -1;
	int sampling = -1;
	int lightType = -1;
    int intensity = -1;

	double fov = PI/4;

	bool surfaceRadiance = true;
	bool mediumRadiance = true;
	bool lightInMedium = true;
	bool cameraInMedium = true;
};

struct Vec {
	double x, y, z;
	Vec(double x0, double y0, double z0) { x = x0; y = y0; z = z0; }
	Vec(double x0) { x = x0; y = x0; z = x0; }
	Vec() { x = 0; y = 0; z = 0; }
	Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x*b, y*b, z*b); }
	Vec operator*(const Vec &b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec operator/(double b) const { return Vec(x / b, y / b, z / b); }
	Vec operator/(const Vec &b) const { return Vec(x / b.x, y / b.y, z / b.z); }
	Vec& norm() { return *this = *this * (1 / std::sqrt(x*x + y * y + z * z)); }
	double length() { return std::sqrt(x*x + y * y + z * z); }
	double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
	Vec operator%(const Vec &b) const { return Vec(y*b.z - z * b.y, z*b.x - x * b.z, x*b.y - y * b.x); }
	double dist(const Vec &b) const { return std::sqrt((x - b.x)*(x - b.x) + (y - b.y)*(y - b.y) + (z - b.z)*(z - b.z)); }
	double average() const { return (x + y + z) / 3.0; }
};

// given v1, set v2 and v3 so they form an orthonormal system
// (we assume v1 is already normalized)
void ons(const Vec& v1, Vec& v2, Vec& v3) {
	if (std::abs(v1.x) > std::abs(v1.y)) {
		// project to the y = 0 plane and construct a normalized orthogonal vector in this plane
		double invLen = 1.f / sqrtf(v1.x * v1.x + v1.z * v1.z);
		v2 = Vec(-v1.z * invLen, 0.0f, v1.x * invLen);
	} else {
		// project to the x = 0 plane and construct a normalized orthogonal vector in this plane
		double invLen = 1.0f / sqrtf(v1.y * v1.y + v1.z * v1.z);
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

// struct to store segments, photon power at the origin of the ray and the distances
struct Segment {
	public:
	Vec origin, direction, end, photonPower;
	double distance;

	Segment(Vec o0, Vec d0, Vec e0, Vec power0, double dist0) :
		origin(o0), direction(d0), end(e0), photonPower(power0), distance(dist0) {}
	Segment() {}
	void changeEnd(Vec e0, double dist0) {
		end = e0; distance = dist0;
	}
	Ray getRay() { return Ray(origin, direction); }
	void setPhotonPower(Vec p) { photonPower = p; }
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

class RestrictedPlane : public Obj {
	public:
	Vec n, min, max;
	double d;
	RestrictedPlane(double d_ = 0, Vec n_ = 0, Vec min_ = 0, Vec max_ = 0) {
		d = d_;
		n = n_;
		min = min_;
		max = max_;
	}
	double intersect(const Ray& ray) const {
		double d0 = n.dot(ray.d);
		if (d0 != 0) {
			double t = -1 * (((n.dot(ray.o)) + d) / d0);

			Vec p;
			if (t > eps) {
				p = ray.o + ray.d*t;
			} else {
				return 0;
			}

			if (n.z != 0) {
				if (min.x < p.x && p.x < max.x &&
					min.y < p.y && p.y < max.y) {
					return t;
				} else {
					return 0;
				}
			} else if (n.y != 0) {
				if (min.x < p.x && p.x < max.x &&
					min.z < p.z && p.z < max.z) {
					return t;
				} else {
					return 0;
				}
			} else if (n.x != 0) {
				if (min.y < p.y && p.y < max.y &&
					min.z < p.z && p.z < max.z) {
					return t;
				} else {
					return 0;
				}
			} else {
				return 0;
			}

			return 0;
		} else return 0;
	}
	Vec normal(const Vec& p0) const { return n; }

	Vec getRandomPoint(double u1, double u2) {
		if (n.z != 0) {
			double distanceX = max.x - min.x;
			double distanceY = max.y - min.y;
			return Vec(min.x + distanceX * u1, min.y + distanceY * u2, -n.z*d);
		} else if (n.y != 0) {
			double distanceX = max.x - min.x;
			double distanceZ = max.z - min.z;
			return Vec(min.x + distanceX * u1, -n.y*d, min.z + distanceZ * u2);
		} else if (n.x != 0) {
			double distanceY = max.y - min.y;
			double distanceZ = max.z - min.z;
			return Vec(-n.x*d, min.y + distanceY * u1, min.z + distanceZ * u2);
		}
		return Vec();
	}
};

class AABox : public Obj {
	public:
	Vec min, max;
	AABox(Vec min_ = 0, Vec max_ = 0) {
		min = min_;
		max = max_;
	}

	//intersection routine from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
	double intersect(const Ray& ray) const {
		double tmin = (min.x - ray.o.x) / ray.d.x;
		double tmax = (max.x - ray.o.x) / ray.d.x;

		if (tmin > tmax) swap(tmin, tmax);

		double tymin = (min.y - ray.o.y) / ray.d.y;
		double tymax = (max.y - ray.o.y) / ray.d.y;

		if (tymin > tymax) swap(tymin, tymax);

		if ((tmin > tymax) || (tymin > tmax))
			return 0;

		if (tymin > tmin)
			tmin = tymin;

		if (tymax < tmax)
			tmax = tymax;

		double tzmin = (min.z - ray.o.z) / ray.d.z;
		double tzmax = (max.z - ray.o.z) / ray.d.z;

		if (tzmin > tzmax) swap(tzmin, tzmax);

		if ((tmin > tzmax) || (tzmin > tmax))
			return 0;

		if (tzmin > tmin)
			tmin = tzmin;

		if (tzmax < tmax)
			tmax = tzmax;

		if (tmin > tmax)
			swap(tmin, tmax);

		return tmin > eps ? tmin : tmax > eps ? tmax : 0;
	}

	Vec normal(const Vec& p0) const {
		if (p0.x >= min.x - eps && p0.x <= min.x + eps) return Vec(-1, 0, 0);
		if (p0.x >= max.x - eps && p0.x <= max.x + eps) return Vec(1, 0, 0);
		if (p0.y >= min.y - eps && p0.y <= min.y + eps) return Vec(0, -1, 0);
		if (p0.y >= max.y - eps && p0.y <= max.y + eps) return Vec(0, 1, 0);
		if (p0.z >= min.z - eps && p0.z <= min.z + eps) return Vec(0, 0, -1);
		if (p0.z >= max.z - eps && p0.z <= max.z + eps) return Vec(0, 0, 1);
	}
};

class Sphere : public Obj {
	public:
	Vec c;
	double r;

	Sphere(double r_ = 0, Vec c_ = 0) { c = c_; r = r_; }
	double intersect(const Ray& ray) const {
		double b = ((ray.o - c) * 2).dot(ray.d);
		double c_ = (ray.o - c).dot((ray.o - c)) - (r*r);
		double disc = b * b - 4 * c_;
		if (disc < 0) return 0;
		else disc = sqrt(disc);
		double sol1 = -b + disc;
		double sol2 = -b - disc;
		return (sol2 > eps) ? sol2 / 2 : ((sol1 > eps) ? sol1 / 2 : 0);
	}

	Vec normal(const Vec& p0) const {
		return (p0 - c).norm();
	}
};

class Intersection {
	public:
	Intersection() { t = inf; object = nullptr; }
	Intersection(double t_, Obj* object_) { t = t_; object = object_; }
	operator bool() { return object != nullptr; }
	double t;
	Obj* object;
};

class Medium {
	public:
	double g, k;
	Vec sigma_a, sigma_s, sigma_t;
	Obj* obj;

	Medium(Vec s_a, Vec s_s, double g0, Obj* o0)
		: sigma_a(s_a), sigma_s(s_s), g(g0), obj(o0) {
		sigma_t = sigma_a + sigma_s;
		k = 1.55*g - 0.55*g*g*g;
	}

	Intersection intersect(const Ray& ray) const {
		Intersection closestIntersection;
		if (!obj) return closestIntersection;
		double t = obj->intersect(ray);
		if (t > eps && t < closestIntersection.t) {
			closestIntersection.t = t;
			closestIntersection.object = obj;
		}
		return closestIntersection;
	}
};

class Pointlight {
	public:
	Vec pos, intensity;
	Pointlight(Vec pos0, Vec intensity0)
		: pos(pos0), intensity(intensity0) {}
};

class Scene {
	vector<Obj*> objects;

	public:
	Medium * medium;
	Pointlight* pointlight;
	double refr_index;

	void add(Obj* object) {
		objects.push_back(object);
	}

	void add(Pointlight* p) {
		pointlight = p;
	}

	void add(Medium* m) {
		medium = m;
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

	void clear() {
		for (int i = 0; i < objects.size(); i++) {
			delete objects.at(i);
		}
	}
};

// Input is the pixel offset, output is the appropriate coordinate
// on the image plane
Vec camcr(const double x, const double y, double fov) {
	double w = width;
	double h = height;
	double fovx = fov;
	double fovy = (h / w) * fovx;
	return Vec(((2 * x - w) / w) * tan(fovx),
			   -((2 * y - h) / h) * tan(fovy),
			   -1.0);
}

double arccos(double x) {
	return acos(x < -1.0 ? -1.0 : x > 1.0 ? 1.0 : x);
}

// rotate a Vector v around another Vector n
Vec rotateVec(Vec n, Vec v, double theta) {
	return v * cos(theta) + n % v*sin(theta);
}

// Uniform sampling on a hemisphere to produce outgoing ray directions.
// courtesy of http://www.rorydriscoll.com/2009/01/07/better-sampling/
Vec sampleHemisphere(double u1, double u2) {
	double r = sqrt(1.0 - u1 * u1);
	double phi = 2 * PI * u2;
	return Vec(cos(phi)*r, sin(phi)*r, u1);
}

Vec sampleSphere(double u1, double u2) {
	const double cos_theta = u1 * 2.0 - 1.0; // remap to -1 .. 1
	const double r = sqrt(1.0 - cos_theta * cos_theta);
	const double phi = 2 * PI * u2;
	return Vec(cos(phi)*r, sin(phi)*r, cos_theta);
}

Vec transmittance(Vec sigma_t, double dist) {
	Vec exp_term = sigma_t * dist*-1;
	return Vec(exp(exp_term.x), exp(exp_term.y), exp(exp_term.z));
}

Vec sampleHenyeyGreenstein(Vec incomingDir, double g, double u1, double u2) {
	double cos_theta;
	if (g == 0) {
		cos_theta = u1 * 2.0 - 1.0;
	} else {
		double g_spu = g * g;
		double squ_term = (1.0 - g_spu) / (1.0 - g + 2 * g* u1);
		cos_theta = (1.0 + g_spu - squ_term * squ_term) / (2.0 * g);
	}
	double sin_theta = sqrt(max(0.0, 1.0 - cos_theta * cos_theta));
	double phi = 2 * PI * u2;
	Vec v1, v2;
	ons(incomingDir, v1, v2);
	return v1 * sin_theta * cos(phi) + v2 * sin_theta * sin(phi) + incomingDir * cos_theta;
}

double pfHenyeyGreenstein(Vec v1, Vec v2, double g) {
	if (g == 0.0) 
		return INV_4PI;
	double cos_theta = v1.dot(v2);
	double g_squ = g * g;
	double temp = 1.0 + g_squ - 2.0*g*cos_theta;
	return INV_4PI * ((1.0 - g_squ) / (temp*sqrt(temp)));
}

//include closestPointOnLine & shortestPointsOnLineSegments functions (just to make this file smaller)
#include "smallpaint_vrl_distanceFunctions.h"

Vec getUniformSample(Segment segment, double &pdf) {
	pdf = 1.0 / segment.distance;
	return segment.origin + segment.direction * segment.distance * RND2;
}

Vec getExponentialSample(Segment segment, double sigma_t, double &pdf) {
	// remap u to account for finite max distance
	double minU = exp(-(sigma_t)* segment.distance);
	double a = RND2 * (1.0 - minU) + minU;

	// sample with pdf proportional to exp(-sig*d)
	double x = -log(a) / (sigma_t);
	pdf = (sigma_t)* a / (1.0 - minU);

	return segment.origin + segment.direction*x;
}

// sampling strategy from Virtual Ray Lights
Vec getLightRaySample(Segment segment, Vec closestPoint, double dist_rayOrigin_closestPoint, double dist_both_rays, double theta_rayDirections, double &pdf) {
	if (theta_rayDirections == 0.0) {
		pdf = 1.0;
		return closestPoint;
	}
	double theta_o = asinh(-dist_rayOrigin_closestPoint / dist_both_rays * sin(theta_rayDirections));
	double theta_e = asinh((segment.distance - dist_rayOrigin_closestPoint) / dist_both_rays * sin(theta_rayDirections));
	double xi = RND2;
	double dist_sample_closestPoint = dist_both_rays * sinh((1 - xi)*theta_o + xi * theta_e) / sin(theta_rayDirections);
	pdf = sin(theta_rayDirections) / ((theta_e - theta_o) * sqrt(dist_both_rays*dist_both_rays + dist_sample_closestPoint * dist_sample_closestPoint * sin(theta_rayDirections)*sin(theta_rayDirections)));
	return  closestPoint + segment.direction*dist_sample_closestPoint;
}

// equi-angular distance sampling strategy from Kulla et al. (Importance sampling of area lights in participating media, 2011)
Vec getEquiAngularSample(Segment segment, Vec closestPoint, double dist_rayOrigin_closestPoint, double dist_interestPoint_closestPoint, double &pdf) {
	double theta_o = atan(-dist_rayOrigin_closestPoint / dist_interestPoint_closestPoint);
	double theta_e = atan((segment.distance - dist_rayOrigin_closestPoint) / dist_interestPoint_closestPoint);
	double xi = RND2;
	double dist_sample_closestPoint = dist_interestPoint_closestPoint * tan((1 - xi)*theta_o + xi * theta_e);
	pdf = dist_interestPoint_closestPoint / ((theta_e - theta_o) * (dist_interestPoint_closestPoint*dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
	return closestPoint + segment.direction*dist_sample_closestPoint;
}

const int M = 10;
struct PiecewiseSegments {
	double k[M - 1], d[M - 1];
	double cdf[M - 1];
};

// compute the slope-intercept form of the current segment as well as the integral
double computeSegment(int i, double theta_s, double theta_e, double fs_s, double fs_e, double &k, double &d, double &cdf) {
	k = (fs_e - fs_s) / (theta_e - theta_s);
	d = fs_e - k * theta_e;
	double integral_end = (k * theta_e * theta_e) / 2.0 + d * theta_e;
	double integral_start = (k * theta_s * theta_s) / 2.0 + d * theta_s;
	cdf = integral_end - integral_start;
	return cdf;
}

// advance sampling proportional to the product of the phase functions and the inverse squared distance
Vec getAngularDomainSample(Segment segment, Vec closestPoint, Vec sample, Vec lightDirection, double g, double dist_rayOrigin_closestPoint, double dist_interestPoint_closestPoint, double &pdf) {
	double theta[M], fs[M];
	PiecewiseSegments segments;

	Vec segmentDirection = segment.direction*-1;

	// construct all variables to find the peak
	Vec a = (segment.end - sample).norm();
	Vec b = (segment.origin - sample).norm();
	Vec c = (a % b).norm();
	// the projection of the light ray on the plane that is spanned by a and b
	Vec e = ((c % lightDirection) % c).norm();

	dist_rayOrigin_closestPoint = segment.distance - dist_rayOrigin_closestPoint;

	// evaluate the angle to the start and the product of the phase functions with this angle
	double	theta_1 = atan(-dist_rayOrigin_closestPoint / dist_interestPoint_closestPoint);
	theta[0] = theta_1;
	fs[0] = pfHenyeyGreenstein(segmentDirection, a, g) * pfHenyeyGreenstein(lightDirection, a, g);

	// evaluate the angle to the end and the product of the phase functions with this angle
	double theta_M = atan((segment.distance - dist_rayOrigin_closestPoint) / dist_interestPoint_closestPoint);
	theta[M - 1] = theta_M;	
	fs[M - 1] = pfHenyeyGreenstein(segmentDirection, b, g) * pfHenyeyGreenstein(lightDirection, b, g);

	double theta_peak;
	bool inRange = false;

	// check if the peak lies on the arc
	if ((a%e).dot((a%b)) >= 0 && (b%e).dot((b%a)) >= 0) {
		// calculate the angle to the peak
		theta_peak = arccos(a.dot(e)) + theta_1;
		inRange = true;
	} else {
		// if the peak is not on the arc, search for the negative peak
		e = e*-1;
		if ((a%e).dot((a%b)) >= 0 && (b%e).dot((b%a)) >= 0) {
			// calculate the angle to the peak
			theta_peak = arccos(a.dot(e)) + theta_1;
			inRange = true;
		}
	}
	int j_peak;
	double cdf = 0.0;

	// case separation for if the peak lies within the arc
	if (inRange) {
		// find the index of the peak
		// this formula form the paper has been adjusted to fit the current indexing
		j_peak = floor((theta_peak - theta_1) / (theta_M - theta_1) * (M - 1) + .5);
		j_peak = j_peak < 2 ? 2 : j_peak;

		// calculate the first sub interval
		for (int j = 2; j <= j_peak; j++) {
			// distribute the angles with a cosine warped uniform spacing
			double theta_j = (theta_peak - theta_1) / 2.0 * (1 - cos((PI*(j - 1)) / (j_peak - 1)));
			theta[j - 1] = theta_j + theta_1;
			// get the direction by rotating the starting direction by the given angle
			Vec dir = rotateVec(c, a, theta_j);
			// calculate the product of the phase functions
			fs[j - 1] = pfHenyeyGreenstein(segmentDirection, dir, g) * pfHenyeyGreenstein(lightDirection, dir, g);
			// compute the slope-intercept form of the current linear function and sum up the intergral/cdf
			cdf += computeSegment(j - 2, theta[j - 2], theta[j - 1], fs[j - 2], fs[j - 1], segments.k[j - 2], segments.d[j - 2], segments.cdf[j - 2]);
		}

		// calculate the second sub interval
		for (int j = j_peak + 1; j <= M - 1; j++) {
			// distribute the angles with a cosine warped uniform spacing
			double theta_j = (theta_M - theta_peak) / 2.0 * (1 - cos((PI*(j - j_peak)) / (M - j_peak)));
			theta[j - 1] = theta_j + theta_peak;
			// get the direction by rotating the starting direction by the given angle
			Vec dir = rotateVec(c, e, theta_j);
			// calculate the product of the phase functions
			fs[j - 1] = pfHenyeyGreenstein(segmentDirection, dir, g) * pfHenyeyGreenstein(lightDirection, dir, g);
			// compute the slope-intercept form of the current linear function and sum up the intergral/cdf
			cdf += computeSegment(j - 2, theta[j - 2], theta[j - 1], fs[j - 2], fs[j - 1], segments.k[j - 2], segments.d[j - 2], segments.cdf[j - 2]);
		}
		// add the last function from the last sampled angle to the end
		cdf += computeSegment(M - 2, theta[M - 2], theta[M - 1], fs[M - 2], fs[M - 1], segments.k[M - 2], segments.d[M - 2], segments.cdf[M - 2]);

	} else {
		// the peak was not in the interval, therefore distribute the angle in just one interval
		for (int j = 2; j <= M - 1; j++) {
			// distribute the angles with a cosine warped uniform spacing
			double theta_j = (theta_M - theta_1) / 2.0 * (1 - cos((PI*(j - 1)) / (M - 1)));
			theta[j - 1] = theta_j + theta_1;
			// get the direction by rotating the starting direction by the given angle
			Vec dir = rotateVec(c, a, theta_j);
			// calculate the product of the phase functions
			fs[j - 1] = pfHenyeyGreenstein(segmentDirection, dir, g) * pfHenyeyGreenstein(lightDirection, dir, g);
			// compute the slope-intercept form of the current linear function and sum up the intergral/cdf
			cdf += computeSegment(j - 2, theta[j - 2], theta[j - 1], fs[j - 2], fs[j - 1], segments.k[j - 2], segments.d[j - 2], segments.cdf[j - 2]);
		}
		// add the last function from the last sampled angle to the end
		cdf += computeSegment(M - 2, theta[M - 2], theta[M - 1], fs[M - 2], fs[M - 1], segments.k[M - 2], segments.d[M - 2], segments.cdf[M - 2]);
	}

	// adjist the xi by the cdf of the piecewise function to be able to find the right function for the evalutaion
	double xi_adj = RND2 * cdf;
	double theta_sample = -1.0;
	double cdf_comb = 0.0;

	// loop through all piecewise segments
	for (int i = 0; i < M-1; i++) {
		// see if the sampled value is in the interval of the cdf for this current segment
		if (cdf_comb < xi_adj && xi_adj <= cdf_comb + segments.cdf[i]) {
			// for a segment that is constant (d = 0) the formula can be simplified
			// this is always the case in isotropic media, but can also occur in anisotropic media
			if (-eps < segments.k[i] && segments.k[i] < eps) {
				// get the sampled angle with the inverse cdf
				theta_sample = (theta[i] * segments.d[i] + xi_adj - cdf_comb) / segments.d[i];
				// calcualte the cdf
				pdf = segments.d[i] / cdf;
			} else {
				// get the sampled angle with the inverse cdf
				double sqrt_term = segments.k[i]*(theta[i] * theta[i] * segments.k[i] + 2.0*xi_adj - 2.0*cdf_comb) + 2.0*theta[i] * segments.d[i]*segments.k[i] + segments.d[i]*segments.d[i];
				theta_sample = (sqrt(sqrt_term) - segments.d[i]) / segments.k[i];
				// calcualte the cdf
				pdf = (segments.k[i]*theta_sample + segments.d[i]) / cdf;
			}
			break;
		}
		cdf_comb += segments.cdf[i];
	}

	// calculate the actual sample
	double dist_sample_closestPoint = dist_interestPoint_closestPoint * tan(theta_sample);
	Vec camSample = closestPoint + segmentDirection *dist_sample_closestPoint;

	// in a case where the pdf is exactly 1, the pdf is equivalent to the equi-angular case
	if (1.0 - eps < pdf && pdf < 1.0 + eps) {
		pdf = dist_interestPoint_closestPoint / ((theta_M - theta_1)*(dist_interestPoint_closestPoint*dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
	} else {
		// the pdf is multiplied by the already adjusted part for the inverse squared distance
		pdf = pdf * (dist_interestPoint_closestPoint / (dist_interestPoint_closestPoint*dist_interestPoint_closestPoint + dist_sample_closestPoint * dist_sample_closestPoint));
	}

	return camSample;
}


// calclualte the radiance from the medium
Vec mediumRadiance(const Scene& scene, std::vector<Segment> lightRays, Segment cameraRay, RenderSettings settings) {
	Vec contribution;
	for (auto iter = lightRays.begin(); iter != lightRays.end(); ++iter) {
		Segment lightRay = *iter;

		// find the points with the shortest distance between the two ray segments
		double dist_both_rays, dist_lightOrigin_closestPointL, dist_cameraOrigin_closestPointC = 0.0;
		Vec closestPointL, closestPointC;
		shortestPointsOnLineSegments(lightRay, cameraRay, dist_both_rays, closestPointL, dist_lightOrigin_closestPointL, closestPointC, dist_cameraOrigin_closestPointC);
		double theta_rayDirections = arccos((lightRay.direction).dot(cameraRay.direction*-1));

		double pdf_light, dist_lightSample_closestPointL, dist_lightOrigin_lightSample;
		Vec lightSample;

		// get the light sample and pdf
		if (settings.sampling == SAMPLE_UNI) {
			lightSample = getUniformSample(lightRay, pdf_light);
		} else if (settings.sampling == SAMPLE_EXP) {
			lightSample = getExponentialSample(lightRay, scene.medium->sigma_t.x, pdf_light);
		} else if (settings.sampling == SAMPLE_ORJ) {
			lightSample = getUniformSample(lightRay, pdf_light);
		} else if (settings.sampling == SAMPLE_SIJ || settings.sampling == SAMPLE_ADJ) {
			lightSample = getLightRaySample(lightRay, closestPointL, dist_lightOrigin_closestPointL, dist_both_rays, theta_rayDirections, pdf_light);
		}

		dist_lightOrigin_lightSample = (lightRay.origin - lightSample).length();

		// find the closest point on the camera ray to the light ray sample
		closestPointOnLine(lightSample, cameraRay, closestPointC, dist_cameraOrigin_closestPointC);
		double dist_lightSample_closestPointC = lightSample.dist(closestPointC);

		// get the sample for the camera ray
		double pdf_camera, dist_cameraSample_closestPoint, dist_cameraOrigin_cameraSample;
		Vec cameraSample;

		// get the camera sample and pdf
		if (settings.sampling == SAMPLE_UNI) {
			cameraSample = getUniformSample(cameraRay, pdf_camera);
		} else if (settings.sampling == SAMPLE_EXP) {
			cameraSample = getExponentialSample(cameraRay, scene.medium->sigma_t.x, pdf_camera);
		} else if (settings.sampling == SAMPLE_ORJ) {
			cameraSample = getEquiAngularSample(cameraRay, closestPointC, dist_cameraOrigin_closestPointC, dist_lightSample_closestPointC, pdf_camera);
		} else if (settings.sampling == SAMPLE_SIJ) {
			cameraSample = getEquiAngularSample(cameraRay, closestPointC, dist_cameraOrigin_closestPointC, dist_lightSample_closestPointC, pdf_camera);
		} else if (settings.sampling == SAMPLE_ADJ) {
			cameraSample = getAngularDomainSample(cameraRay, closestPointC, lightSample, lightRay.direction, scene.medium->g, dist_cameraOrigin_closestPointC, dist_lightSample_closestPointC, pdf_camera);
		}

		dist_cameraOrigin_cameraSample = (cameraRay.origin - cameraSample).length();

		// calculate distance and direction between samples
		Vec dir_samples = (lightSample - cameraSample).norm();
		double dist_samples = cameraSample.dist(lightSample);

		// calculating visibility
		Intersection intersection = scene.intersect(Ray(cameraSample, dir_samples));
		if (intersection) {
			if (intersection.t < dist_samples - eps) {
				continue;
			}
		}
		// phase functions with the angle as input (Henyey-Greenstein phase function)
		double fs_u = pfHenyeyGreenstein(cameraRay.direction, dir_samples, scene.medium->g);
		double fs_v = pfHenyeyGreenstein(lightRay.direction*-1, dir_samples, scene.medium->g);

		// transmittance functions
		Vec trans = transmittance(scene.medium->sigma_t, dist_cameraOrigin_cameraSample + dist_lightOrigin_lightSample + dist_samples);

		// inverse squared distance term
		double geom = 1.0 / (dist_samples*dist_samples);

		// putting the integrand of the formula together
		Vec integrand = lightRay.photonPower *scene.medium->sigma_s*scene.medium->sigma_s*trans*fs_u*fs_v*geom;
		contribution = contribution + (integrand / (pdf_light * pdf_camera)) / settings.vrlps;
	}

	return contribution;
}

// calculate the surface radiance
Vec surfaceRadiance(const Scene& scene, std::vector<Segment> lightRays, Segment cameraRay, Vec surfaceColor, Vec surfaceNormal, RenderSettings settings) {
	Vec contribution;
	for (auto iter = lightRays.begin(); iter != lightRays.end(); ++iter) {
		Segment lightRay = *iter;

		// find the closest point on the light ray to the camera ray.e
		double dist_cameraEnd_closestPointL, dist_lightOrigin_closestPointL;
		Vec closestPointL;
		closestPointOnLine(cameraRay.end, lightRay, closestPointL, dist_lightOrigin_closestPointL);
		dist_cameraEnd_closestPointL = cameraRay.end.dist(closestPointL);

		double pdf_light, dist_lightSample_closestPoint, dist_lightOrigin_lightSample;
		Vec lightSample;

		// get the light sample and pdf
		if (settings.sampling == SAMPLE_UNI) {
			lightSample = getUniformSample(lightRay, pdf_light);
		} else if (settings.sampling == SAMPLE_EXP) {
			lightSample = getExponentialSample(lightRay, scene.medium->sigma_t.x, pdf_light);
		} else if (settings.sampling == SAMPLE_ORJ) {
			lightSample = getEquiAngularSample(lightRay, closestPointL, dist_lightOrigin_closestPointL, dist_cameraEnd_closestPointL, pdf_light);
		} else if (settings.sampling == SAMPLE_SIJ) {
			lightSample = getEquiAngularSample(lightRay, closestPointL, dist_lightOrigin_closestPointL, dist_cameraEnd_closestPointL, pdf_light);
		} else if (settings.sampling == SAMPLE_ADJ) {
			lightSample = getAngularDomainSample(lightRay, closestPointL, cameraRay.end, cameraRay.direction, scene.medium->g, dist_lightOrigin_closestPointL, dist_cameraEnd_closestPointL, pdf_light);
		}

		dist_lightOrigin_lightSample = (lightRay.origin - lightSample).length();

		// calculate distance and direction between samples
		Vec dir_sample = (lightSample - cameraRay.end).norm();
		double dist_sample = lightSample.dist(cameraRay.end);

		// calculating visibility
		Intersection intersection = scene.intersect(Ray(cameraRay.end, dir_sample));
		if (intersection) {
			if (intersection.t < dist_sample - eps) {
				continue;
			}
		}

		// phase functions with the angle as input (Henyey-Greenstein phase function)
		double fs_v = pfHenyeyGreenstein(lightRay.direction*-1, dir_sample, scene.medium->g);

		// cosine-weight the BRDF
		Vec cosWeightedSurfaceColor = surfaceColor * abs(surfaceNormal.dot(dir_sample*-1));

		// transmittance functions
		Vec trans = transmittance(scene.medium->sigma_t, dist_sample + dist_lightOrigin_lightSample);

		// inverse squared distance term
		double geom = 1.0 / (dist_sample*dist_sample);

		// putting the integrand of the formula together
		Vec integrand = lightRay.photonPower*cosWeightedSurfaceColor*scene.medium->sigma_s*fs_v*trans*geom;

		contribution = contribution + (integrand / pdf_light) / settings.vrlps;
	}
	return contribution;
}

void trace(Segment segment, const Scene& scene, std::vector<Segment>& lightRays, bool insideMedium, int depth, Vec& clr, RenderSettings settings, int pass) {

	double freeFlightDistance = 0.0;
	if (pass == LIGHTPASS) {
		if (depth >= settings.bounces)
			return;

		double freeFlightDistance = -log(RND2) / scene.medium->sigma_t.average();

		if (!insideMedium)
			freeFlightDistance = inf;

	}

	double rrFactor = 1.0;
	if (pass == RENDERPASS) {
		// Russian roulette: starting at depth 3, each recursive step will stop with a probability of 0.1
		if (depth >= 3) {
			const double rrStopProbability = 0.1;
			if (RND2 <= rrStopProbability) {
				return;
			}
			rrFactor = 1.0 / (1.0 - rrStopProbability);
		}
	}

	Intersection intersection = scene.intersect(segment.getRay());
	Intersection intersectionMedium = scene.medium->intersect(segment.getRay());

	bool mediumTransition = false;
	bool scatteringEvent = false;

	if (freeFlightDistance < intersection.t && freeFlightDistance < intersectionMedium.t)
		scatteringEvent = true;

	// a medium transition can only occur at an intersection with the medium where there is no intersection before this event
	if (intersectionMedium) {
		if (!intersection) {
			mediumTransition = true;
		} else if (intersectionMedium.t < intersection.t) {
			mediumTransition = true;
		}
	}

	Segment original = segment;

	// no intersection with geometry, no medium transition and outside the medium
	if (!intersection && !mediumTransition && !insideMedium) return;

	// no intersection with geometry, no medium transition and inside the medium
	if (!intersection && !mediumTransition && insideMedium) {
		segment.changeEnd(segment.origin + segment.direction * inf, inf);
		if (pass == LIGHTPASS) {
			lightRays.push_back(segment);
		}
		if (pass == RENDERPASS) {
			if (settings.mediumRadiance) clr = clr + mediumRadiance(scene, lightRays, segment, settings) * rrFactor;
			return;
		}
	}

	// intersection, but no medium transition
	if (intersection && !mediumTransition) {
		// Travel the ray to the hit point where the closest object lies and compute the surface normal there.
		Vec hp = segment.origin + segment.direction * intersection.t;
		Vec N = intersection.object->normal(hp);
		Vec surfaceClr = intersection.object->cl;
		int type = intersection.object->type;


		if (pass == LIGHTPASS) {
			segment.changeEnd(hp, intersection.t);
			if (insideMedium) lightRays.push_back(segment);
			if (depth + 1 >= settings.bounces) return;
		}

		if (pass == RENDERPASS) {
			// calculate VRL contribution
			if (insideMedium) {
				segment.changeEnd(hp, intersection.t);
				if (settings.mediumRadiance) clr = clr + mediumRadiance(scene, lightRays, segment, settings) * rrFactor;
			} else {
				segment.changeEnd(hp, 0);
			}
		}

		// Diffuse BRDF
		if (type == DIFFUSE) {
			if (pass == LIGHTPASS) {
				// create new ray direction
				Vec rotX, rotY;
				ons(N, rotX, rotY);
				Vec sampledDir = sampleHemisphere(RND2, RND2);
				Vec rotatedDir;
				rotatedDir.x = Vec(rotX.x, rotY.x, N.x).dot(sampledDir);
				rotatedDir.y = Vec(rotX.y, rotY.y, N.y).dot(sampledDir);
				rotatedDir.z = Vec(rotX.z, rotY.z, N.z).dot(sampledDir);
				Vec direction = rotatedDir;	// already normalized

				// calculate the diffuse reflection probability
				double diff_reflect_prob = max(surfaceClr.x, max(surfaceClr.y, surfaceClr.z));
				// calculate the new photon photonPower
				Vec photonPower;
				photonPower.x = segment.photonPower.x * (surfaceClr.x / diff_reflect_prob);
				photonPower.y = segment.photonPower.y * (surfaceClr.y / diff_reflect_prob);
				photonPower.z = segment.photonPower.z * (surfaceClr.z / diff_reflect_prob);
				segment.setPhotonPower(photonPower);

				// multiply by transmittance
				if (insideMedium) {
					segment.setPhotonPower(segment.photonPower * transmittance(scene.medium->sigma_t, segment.distance));
				}
				// trace new light ray

				segment = Segment(segment.end, direction, Vec(), segment.photonPower, 0);
				trace(segment, scene, lightRays, insideMedium, depth + 1, clr, settings, LIGHTPASS);
			}

			if (pass == RENDERPASS) {
				// calculate surface radiance and multiply by transmittance and russion roulette factor
				if (settings.surfaceRadiance) clr = clr + surfaceRadiance(scene, lightRays, segment, surfaceClr / 255.0, N, settings) * rrFactor;
			}
		}

		// Specular BRDF - this is a singularity in the rendering equation that follows
		// delta distribution, therefore we handle this case explicitly - one incoming
		// direction -> one outgoing direction, that is, the perfect reflection direction.
		if (type == SPECULAR) {
			Vec direction = (segment.direction - N * (segment.direction.dot(N) * 2)).norm();

			if (pass == LIGHTPASS) {
				// multiply by transmittance
				if (insideMedium) {
					segment.setPhotonPower(segment.photonPower * transmittance(scene.medium->sigma_t, segment.distance));
				}
				// create new light ray
				segment = Segment(segment.end, direction, Vec(), segment.photonPower, 0);
				trace(segment, scene, lightRays, insideMedium, depth + 1, clr, settings, LIGHTPASS);
			}

			if (pass == RENDERPASS) {
				segment = Segment(segment.end, direction, Vec(), Vec(), 0);
				Vec temp;
				trace(segment, scene, lightRays, insideMedium, depth + 1, temp, settings, RENDERPASS);
				clr = clr + temp;
			}
		}


		// Glass/refractive BRDF - we use the vector version of Snell's law and Fresnel's law
		// to compute the outgoing reflection and refraction directions and probability weights.
		if (type == REFRACTIVE) {
			Vec direction;
			double n = scene.refr_index;
			double R0 = (1.0 - n) / (1.0 + n);
			R0 = R0 * R0;
			if (N.dot(segment.direction) > 0) { // we're inside the medium
				N = N * -1;
				n = 1 / n;
			}
			n = 1 / n;
			double cost1 = (N.dot(segment.direction))*-1; // cosine of theta_1
			double cost2 = 1.0 - n * n*(1.0 - cost1 * cost1); // cosine of theta_2
			double Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0); // Schlick-approximation
			if (cost2 > 0 && RND2 > Rprob) { // refraction direction
				direction = ((segment.direction*n) + (N*(n*cost1 - sqrt(cost2)))).norm();
			} else { // reflection direction
				direction = (segment.direction + N * (cost1 * 2)).norm();
			}

			if (pass == LIGHTPASS) {
				// multiply by transmittance
				if (insideMedium) {
					segment.setPhotonPower(segment.photonPower * transmittance(scene.medium->sigma_t, segment.distance));
				}
				// create new light ray
				segment = Segment(segment.end, direction, Vec(), segment.photonPower, 0);
				trace(segment, scene, lightRays, insideMedium, depth + 1, clr, settings, LIGHTPASS);
			}

			if (pass == RENDERPASS) {
				segment = Segment(segment.end, direction, Vec(), Vec(), 0);
				Vec temp;
				trace(segment, scene, lightRays, insideMedium, depth + 1, temp, settings, RENDERPASS);
				clr = clr + temp;
			}
		}
	} else if (mediumTransition) {
		Vec mediumTransitionPoint = segment.origin + segment.direction * intersectionMedium.t;


		if (pass == LIGHTPASS) {
			segment.changeEnd(mediumTransitionPoint, intersectionMedium.t);

			if (insideMedium) {
				lightRays.push_back(segment);
				segment.setPhotonPower(segment.photonPower * transmittance(scene.medium->sigma_t, segment.distance));
			}

			segment = Segment(segment.end, segment.direction, Vec(), segment.photonPower, 0);
			trace(segment, scene, lightRays, !insideMedium, depth, clr, settings, LIGHTPASS);
		}

		if (pass == RENDERPASS) {
			if (insideMedium) {
				segment.changeEnd(mediumTransitionPoint, intersectionMedium.t);
				if (settings.mediumRadiance) clr = clr + mediumRadiance(scene, lightRays, segment, settings) * rrFactor;
			} else {
				segment.changeEnd(mediumTransitionPoint, 0);
			}
			segment = Segment(segment.end, segment.direction, Vec(), Vec(), 0);
			Vec temp;
			trace(segment, scene, lightRays, !insideMedium, depth, temp, settings, RENDERPASS);
			clr = clr + temp;
		}
	}

	if (pass == LIGHTPASS && scatteringEvent) {
		if (depth >= settings.bounces - 1)
			return;

		segment = original;

		Vec scatteringPoint = segment.origin + segment.direction * freeFlightDistance;
		segment.changeEnd(scatteringPoint, freeFlightDistance);

		double scattering_prob = (scene.medium->sigma_s / scene.medium->sigma_t).average();
		if (RND2 > scattering_prob) return;

		segment.setPhotonPower(segment.photonPower * transmittance(scene.medium->sigma_t, segment.distance));

		// create new ray direction
		Vec direction = sampleHenyeyGreenstein(segment.direction, scene.medium->g, RND2, RND2);

		// create new light ray
		segment = Segment(segment.end, direction, Vec(), segment.photonPower, 0);
		trace(segment, scene, lightRays, insideMedium, depth + 1, clr, settings, LIGHTPASS);
	}
}

void clear(Scene& scene, Vec **pix) {
	scene.clear();
	for (int i = 0; i < width; i++) {
		delete[] pix[i];
	}
	delete[] pix;
}

int color(double clr) {
	return min((int)(clr * 255.0), 255);
}

#include "smallpaint_vrl_tests.h"

void saveImage(Vec **pix, std::string name, int spp) {
	FILE *f = fopen(std::string(name + ".ppm").c_str(), "w");
	fprintf(f, "P3\n%d %d\n%d\n ", width, height, 255);
	for (int row = 0; row < height; row++) {
		for (int col = 0; col < width; col++) {

			fprintf(f, "%d %d %d ", color(pix[col][row].x / (double)spp), color(pix[col][row].y / (double)spp), color(pix[col][row].z / (double)spp));
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

bool progression = false;
bool timecomp = false;
string testname = "";

void render(int id, int size, int spp, int vrlps, int sceneID, int bounces, double sigma_a, double sigma_s, double g, int sampling, double refr_index, int intensity, bool mediumRadiance, bool surfaceRadiance) {
    RenderSettings settings = RenderSettings();
    settings.spp = spp;
    settings.vrlps = vrlps;
    settings.scene = sceneID;
    settings.bounces = bounces;
    settings.sigma_a = sigma_a;
    settings.sigma_s = sigma_s;
    settings.g = g;
    settings.sampling = sampling;
    settings.intensity = intensity;
    settings.mediumRadiance = mediumRadiance;
    settings.surfaceRadiance = surfaceRadiance;



    Scene scene;
	loadScene(scene, settings);

	// create an area light with four planes if needed
	PlaneAreaLights areaLights;
	if (settings.lightType == LIGHT_DIR) {
		areaLights.add(new RestrictedPlane(2.5 - eps, Vec(0, -1, 0), Vec(-.65, 0, -3.5), Vec(-.495, 0, -3.18)));
		areaLights.add(new RestrictedPlane(2.5 - eps, Vec(0, -1, 0), Vec(-.425, 0, -3.5), Vec(-.27, 0, -3.18)));
		areaLights.add(new RestrictedPlane(2.5 - eps, Vec(0, -1, 0), Vec(-.2, 0, -3.5), Vec(-.045, 0, -3.18)));
		areaLights.add(new RestrictedPlane(2.5 - eps, Vec(0, -1, 0), Vec(.025, 0, -3.5), Vec(.18, 0, -3.18)));
	}

    scene.refr_index = refr_index;

    width = size;
    height = size;

	Vec **pix = new Vec*[width];
	for (int i = 0; i < width; i++) {
		pix[i] = new Vec[height];
	}

	// generate equally spaced direction
	std:vector<Vec> directions;
	fibonacciSamples(settings.vrlps, directions);

    bool running = true;
	int vrls = 0;
	for (int i = 0; i < settings.spp; i++) {

		std::vector<Segment> lightRays;
		// for the number of vrl per sample, cast a light ray in and follow it through the scene
		for (int j = 0; j < settings.vrlps; j++) {
			Vec direction;
			if (settings.lightType == LIGHT_DIR) {
				scene.pointlight->pos = areaLights.getRandomPoint(RND2, RND2, RND2);
				direction = areaLights.getDirection();
			} else if (settings.lightType == LIGHT_BIDIR) {
				direction = directions.at(j);
			} else {
				direction = sampleSphere(RND2, RND2).norm();
			}

			// 
			Segment segment = Segment(scene.pointlight->pos, direction, Vec(), scene.pointlight->intensity, 0);
			Vec temp;
			trace(segment, scene, lightRays, settings.lightInMedium, 0, temp, settings, LIGHTPASS);
		}
		vrls = vrls + lightRays.size();

        #pragma omp parallel for schedule(dynamic)
		for (int col = 0; col < width; col++) {
			for (int row = 0; row < height; row++) {
                if (running) {
                    Vec color;
                    Ray ray;
                    ray.o = (Vec(0, 0, 0)); // rays start out from here
                    Vec cam = camcr(col, row, settings.fov); // construct image plane coordinates
                    cam.x = cam.x + RND / 700; // anti-aliasing for free
                    cam.y = cam.y + RND / 700;
                    ray.d = (cam - ray.o).norm(); // point from the origin to the camera plane
                    Segment segment = Segment(ray.o, ray.d, Vec(), Vec(), 0);

                    trace(segment, scene, lightRays, settings.cameraInMedium, 0, color, settings, RENDERPASS);
                    pix[col][row] = pix[col][row] + color;
                    if (!smallpaint::isRunning(id)) running = false;
                }
			}
		}
        if (!running) return;
        imageOutput(pix, i);
	}


    /*string str_g = to_string(abs(settings.g)).substr(0, 4);
	str_g.erase(remove(str_g.begin(), str_g.end(), '.'), str_g.end());
	std::string name = "vrl_spp_" + to_string(settings.spp) + "_vrlps_" + to_string(settings.vrlps) + "_b_" + to_string(settings.bounces) + (settings.g > 0 ? "_g_+" : settings.g == 0 ? "_g_" : "_g_-") + str_g + "_t_" + to_string(settings.sampling);
	saveImage(pix, name, settings.spp);
	
	clear(scene, pix);
    return 0;*/
}
}
