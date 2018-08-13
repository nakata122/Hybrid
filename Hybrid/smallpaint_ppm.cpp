// smallpaint by karoly zsolnai - zsolnai@cg.tuwien.ac.at
//
// This program is used as an educational learning tool on the Rendering
// course at TU Wien. Course webpage:
// http://cg.tuwien.ac.at/courses/Rendering/
//
// This renderer is a different version of this work
// http://www.ci.i.u-tokyo.ac.jp/~hachisuka/smallppm.cpp
// by Hachisuka. It was modified and implemented into
// Smallpaint by Michael Oppitz - e1227129@student.tuwien.ac.at

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctime>
#include <string>
#include <unordered_map>
#include <random>
#include <cstdint>
#include <algorithm>
#include <sstream>
#include <queue>
#include "../src_gui/main.h"

namespace smallpaint_ppm {

// Helpers for random number generation
std::mt19937 mersenneTwister;
std::uniform_real_distribution<double> uniform;

#define RND (2.0*uniform(mersenneTwister)-1.0)
#define RND2 (uniform(mersenneTwister))

#define ALPHA ((double)0.7)
#define PI ((double)3.14159265358979)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int width, height;
const double inf = 1e9;
const double eps = 1e-6;
using namespace std;
typedef unordered_map<string, double> pl;
unsigned int hash_num;
float hash_size;

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
	double max() const { return x > y&&x > z ? x : y > z ? y : z; }
	Vec minValues(const Vec &b) const { return Vec(MIN(x, b.x), MIN(y, b.y), MIN(z, b.z)); }
	Vec maxValues(const Vec &b) const { return Vec(MAX(x, b.x), MAX(y, b.y), MAX(z, b.z)); }
	Vec operator%(const Vec &b) const { return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }
    Vec absolute() const { return Vec(abs(x), abs(y), abs(z)); }
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
	double emission = 0;
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
		double b = ((ray.o - c) * 2).dot(ray.d);
		double c_ = (ray.o - c).dot((ray.o - c)) - (r*r);
		double disc = b*b - 4 * c_;
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

class Scene {
	vector<Obj*> objects;

public:
	void add(Obj* object) {
		objects.push_back(object);
	}

	void remove(Obj* object) {
		objects.erase(std::remove(objects.begin(), objects.end(), object), objects.end());
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

struct AABB {
	Vec min, max;
	inline void fit(const Vec &p) {
		min = p.minValues(min);
		max = p.maxValues(max);
	}
	inline void reset() {
		min = Vec(inf, inf, inf);
		max = Vec(-inf, -inf, -inf);
	}
};

struct HitPoint {
	Vec baseColor, position, normal, flux;
	double radius;
	unsigned int photonCount;
	int x, y;
	bool light = false;
};

struct List {
	HitPoint *id;
	List *next;
};

List* ListAdd(HitPoint *i, List* h) {
	List* list = new List;
	list->id = i;
	list->next = h;
	return list;
}

List **hash_grid;
List *hitpoints = NULL;
std::vector<HitPoint*> nohit;
std::vector<HitPoint*> veryhit;
AABB hitpointsBB;

namespace hashNS {
inline unsigned int hash(const int ix, const int iy, const int iz) {
	return (unsigned int)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % hash_num;
}
}

void build_hash_grid(const int w, const int h) {
	// add all HitPoints to BB
	hitpointsBB.reset();
	List *lst = hitpoints;
	while (lst != NULL) {
		HitPoint *hp = lst->id;
        lst = lst->next;
        nohit.push_back(hp);
        veryhit.push_back(hp);
		hitpointsBB.fit(hp->position);
	}

	// create initial radius
	Vec ssize = hitpointsBB.max - hitpointsBB.min;
    float irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;

	// add all HitPoints + radius to BB
	hitpointsBB.reset();
	lst = hitpoints;
	while (lst != NULL) {
		HitPoint *hp = lst->id;
		lst = lst->next;
		hp->radius = irad*irad;
		hp->photonCount = 0;
		hp->flux = Vec();
		hash_num++;
		hitpointsBB.fit(hp->position - irad);
		hitpointsBB.fit(hp->position + irad);
	}

	// make each grid cell two times larger than the initial radius
    hash_size = 1.0 / (irad*2.0);

	// build the hash table
	hash_grid = new List*[hash_num];
	for (unsigned int i = 0; i < hash_num; i++) hash_grid[i] = NULL;
	lst = hitpoints;
	while (lst != NULL) {
		HitPoint *hp = lst->id;
		lst = lst->next;

		Vec BMin = ((hp->position - irad) - hitpointsBB.min) * hash_size;
		Vec BMax = ((hp->position + irad) - hitpointsBB.min) * hash_size;

		int xAdd = 0, yAdd = 0, zAdd = 0;
		if (BMin.x == BMax.x) xAdd = 1;
		if (BMin.y == BMax.y) yAdd = 1;
		if (BMin.z == BMax.z) zAdd = 1;

		for (int iz = abs(int(BMin.z)) - zAdd; iz <= abs(int(BMax.z)) + zAdd; iz++) {
			for (int iy = abs(int(BMin.y)) - yAdd; iy <= abs(int(BMax.y)) + yAdd; iy++) {
				for (int ix = abs(int(BMin.x)) - xAdd; ix <= abs(int(BMax.x)) + xAdd; ix++) {
					int hv = hashNS::hash(ix, iy, iz);
					hash_grid[hv] = ListAdd(hp, hash_grid[hv]);
				}
			}
		}
	}
}

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
	const double r = sqrt(u1*u1);
	const double phi = 2 * PI * u2;
	return Vec(cos(phi)*r, sin(phi)*r, sqrt(1 - u1*u1));
}

//returns a random point on the surface
Ray sphere(Obj* light, double u1, double u2) {
	Sphere* s = dynamic_cast<Sphere*>(light);
	double theta = 2 * PI * u1;
	double phi = acos(1 - 2 * u2);
	double x = sin(phi) * cos(theta);
	double y = sin(phi) * sin(theta);
	double z = cos(phi);
	Vec d = Vec(x, y, z);
	return Ray(s->c + d * s->r * (1.1), Vec(x, y, z));
}

void clear(Scene& scene, Vec **pix) {
	List *lst, *next;

	for (unsigned int i = 0; i < hash_num; i++) {
		lst = hash_grid[i];
		while (lst != NULL) {
			next = lst->next;
            delete lst;
			lst = next;
		}
	}
    delete[] hash_grid;

	lst = hitpoints;
	while (lst != NULL) {
		HitPoint* hp = lst->id;
        delete hp;
		next = lst->next;
        delete lst;
		lst = next;
	}

	hitpoints = NULL;
	hitpointsBB.reset();

	scene.clear();
	for (int i = 0; i < width; i++) {
        delete[] pix[i];
	}
    delete[] pix;
}

void trace(Ray &ray, const Scene& scene, int depth, bool hitpointPass, Vec& flux, Vec& baseCl, int pixel_x, int pixel_y, pl& params) {
	// Russian roulette: starting at depth 5, each recursive step will stop with a probability of 0.1
    if (depth >= 3) {
		const double rrStopProbability = 0.1;
		if (RND2 <= rrStopProbability) {
			return;
		}
	}

	Intersection intersection = scene.intersect(ray);
	if (!intersection) return;

	// Travel the ray to the hit point where the closest object lies and compute the surface normal there.
	Vec hp = ray.o + ray.d * intersection.t;
	Vec N = intersection.object->normal(hp);
    ray.o = hp;



    if (intersection.object->type == 1) {

        if (hitpointPass) {
            // eye ray
            HitPoint* hitpoint = new HitPoint;
            if (intersection.object->emission > 0) {
                hitpoint->light = true;
            }
            hitpoint->baseColor = intersection.object->cl.mult(baseCl);
            hitpoint->position = hp;
            hitpoint->normal = N;
            hitpoint->x = pixel_x;
            hitpoint->y = pixel_y;
            hitpoints = ListAdd(hitpoint, hitpoints);
            flux = intersection.object->cl;
        } else {
            // photon ray
            Vec hashLocation = ((hp - hitpointsBB.min) * hash_size).absolute();
        #pragma omp critical
            {
                List* hitpoints = hash_grid[hashNS::hash(hashLocation.x, hashLocation.y, hashLocation.z)];
                while (hitpoints != NULL) {
                    HitPoint *hitpoint = hitpoints->id;
                    hitpoints = hitpoints->next;

                    if (!hitpoint->light) {

                        Vec diff = hitpoint->position - hp;

                        // check normals to be closer than 90 degree (avoids some edge brightning)
                        if ((hitpoint->normal.dot(N) > 1e-3) && (diff.dot(diff) <= hitpoint->radius)) {
                            float alphaCorrectedPhotonCount = hitpoint->photonCount*ALPHA;
                            double radiusCorrection = (alphaCorrectedPhotonCount + ALPHA) / (alphaCorrectedPhotonCount + 1.0);
                            hitpoint->radius = hitpoint->radius*radiusCorrection;
                            hitpoint->photonCount++;
                            float dist = 1 - (intersection.t / (hitpointsBB.max.max() - hitpointsBB.min.max()));
                            Vec newFlux = hitpoint->baseColor.mult(flux)*(1. / PI);
                            hitpoint->flux = (hitpoint->flux + newFlux.mult(Vec(dist, dist, dist)))*radiusCorrection;

                            if (intersection.object->emission > 0) hitpoint->flux = hitpoint->flux + intersection.object->cl * intersection.object->emission *(PI*4.0);
                        }
                    }
                }
            }

            if (RND2 < intersection.object->cl.max()) {
                HitPoint *hip = nohit[RND2*nohit.size()];
                if(acos((hip->position - ray.o).norm().dot(N)) < 90){
                    ray.d = (hip->position - ray.o).norm();
                }else{
                    Vec rotX, rotY;
                    ons(N, rotX, rotY);
                    Vec sampledDir = hemisphere(RND2, RND2);
                    Vec rotatedDir;
                    rotatedDir.x = Vec(rotX.x, rotY.x, N.x).dot(sampledDir);
                    rotatedDir.y = Vec(rotX.y, rotY.y, N.y).dot(sampledDir);
                    rotatedDir.z = Vec(rotX.z, rotY.z, N.z).dot(sampledDir);
                    ray.d = rotatedDir;	// already normalized
                }
                Vec newFlux = intersection.object->cl.mult(flux)*(1. / intersection.object->cl.max());
                trace(ray, scene, depth + 1, hitpointPass, newFlux, baseCl, pixel_x, pixel_y, params);
            }
        }
    } else if (intersection.object->type == 2) {
        double cost = ray.d.dot(N);
        ray.d = (ray.d - N*(cost * 2)).norm();
        Vec newFlux = intersection.object->cl.mult(flux);
        Vec newbaseCl = intersection.object->cl.mult(baseCl);
        trace(ray, scene, depth + 1, hitpointPass, newFlux, newbaseCl, pixel_x, pixel_y, params);
    } else {
        double n = params["refr_index"];
        double R0 = (1.0 - n) / (1.0 + n);
        R0 = R0*R0;
        if (N.dot(ray.d) > 0) { // we're inside the medium
            N = N*-1;
            n = 1 / n;
        }
        n = 1 / n;
        double cost1 = (N.dot(ray.d))*-1; // cosine of theta_1
        double cost2 = 1.0 - n*n*(1.0 - cost1*cost1); // cosine of theta_2
        double Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0); // Schlick-approximation
        if (cost2 > 0 && RND2 > Rprob) { // refraction direction
            ray.d = ((ray.d*n) + (N*(n*cost1 - sqrt(cost2)))).norm();
        } else { // reflection direction
            ray.d = (ray.d + N*(cost1 * 2)).norm();
        }
        Vec newFlux = intersection.object->cl.mult(flux);
        Vec newbaseCl = intersection.object->cl.mult(baseCl);
        trace(ray, scene, depth + 1, hitpointPass, newFlux, newbaseCl, pixel_x, pixel_y, params);
    }
}

void render(int id, int size, int photons, double refr_index) {
	pl params;

	Scene scene;
	auto add = [&scene](Obj* s, Vec cl, double emission, int type) {
		s->setMat(cl, emission, type);
		scene.add(s);
	};

	// Radius, position, color, emission, type (1=diff, 2=spec, 3=refr) for spheres
    add(new Sphere(1.05, Vec(-0.75, -1.45, -4.4)), Vec(1, 1, 1), 0, 2); // Middle sphere
    add(new Sphere(0.5, Vec(2.0, -2.05, -3.7)), Vec(1, 1, 1), 0, 3); // Right sphere
    add(new Sphere(0.6, Vec(-1.75, -1.95, -3.1)), Vec(.01, .01, .5), 0, 1); // Left sphere
	// Position, normal, color, emission, type for planes
    add(new Plane(2.5, Vec(0, 1, 0)), Vec(.1, .1, .1), 0, 1); // Bottom plane
    add(new Plane(5.5, Vec(0, 0, 1)), Vec(.1, .1, .1), 0, 1); // Back plane
    add(new Plane(2.75, Vec(1, 0, 0)), Vec(.5, .01, .01), 0, 1); // Left plane
    add(new Plane(2.75, Vec(-1, 0, 0)), Vec(.01, .5, .01), 0, 1); // Right plane
    add(new Plane(3.0, Vec(0, -1, 0)), Vec(.1, .1, .1), 0, 1); // Ceiling plane
    add(new Plane(0.5, Vec(0, 0, -1)), Vec(.1, .1, .1), 0, 1); // Front plane

	Obj* light = new Sphere(0.5, Vec(0, 1.9, -3));
    light->setMat(Vec(1, 1, 1), 50, 1);
	scene.add(light);

	params["refr_index"] = refr_index;

	width = size;
	height = size;

	// trace eye rays and store measurement points
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Ray ray;
            Vec cam = camcr(x, y); // construct image plane coordinates
            ray.d = (cam - ray.o).norm(); // point from the origin to the camera plane
            Vec temp1 = Vec(1, 1, 1);
            Vec temp2 = Vec(1, 1, 1);
            trace(ray, scene, 0, true, temp1, temp2, x, y, params);
        }
    }

	// build the hash table over the measurement points
	build_hash_grid(width, height);

	Vec **pix = new Vec*[width];
	for (int i = 0; i < width; i++) {
		pix[i] = new Vec[height];
	}

	bool running = true;

	Ray ray;
    Vec flux = light->cl*light->emission*(PI*4.0);

	int updateF = 150;
	if (size >= 500) updateF = 50;
	if (size >= 1000) updateF = 25;

    for (int i = 1; i <= photons / 1000; i++) {
    #pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < 1000; j++) {
			if (running) {
                ray = sphere(light, RND2, RND2);

                if(RND2 >= 0.5){
                    HitPoint *hip = veryhit[RND2*veryhit.size()];
                    ray.d = (hip->position - ray.o).norm();
                }else {
                    Vec rotX, rotY;
                    ons(ray.o, rotX, rotY);
                    Vec sampledDir = hemisphere(RND2, RND2);
                    Vec rotatedDir;
                    rotatedDir.x = Vec(rotX.x, rotY.x, ray.d.x).dot(sampledDir);
                    rotatedDir.y = Vec(rotX.y, rotY.y, ray.d.y).dot(sampledDir);
                    rotatedDir.z = Vec(rotX.z, rotY.z, ray.d.z).dot(sampledDir);
                    ray.d = rotatedDir;	// already normalizedHitPoint
                }
				Vec temp = Vec(1, 1, 1);
                trace(ray, scene, 0, false, flux, temp, -1, -1, params);
				if (!smallpaint::isRunning(id)) running = false;
			}
		}

		// density estimation
        List* lst = hitpoints;

        nohit.resize(2);
        veryhit.resize(2);
        while (lst != NULL) {
            HitPoint* hp = lst->id;
            Vec fl = hp->flux*(1.0 / (PI*hp->radius* i * 1000.0));
            //Every point that needs light
            if(fl.x + fl.y + fl.z <= 0.2 || hp->photonCount <= 1)
                nohit.push_back(hp);
            if(hp->photonCount < 1)
                veryhit.push_back(hp);
            lst = lst->next;
            if (hp->light) {
                pix[hp->x][hp->y] = hp->baseColor * 255;
            } else {
                pix[hp->x][hp->y] = fl;
            }
        }

		if (!running) {
			clear(scene, pix);
			return;
		}

		if (i%updateF == 0) imageOutput(pix, i * 1000 - 1);
	}
	if ((photons / 1000) % updateF != 0) imageOutput(pix, photons - 1);
	clear(scene, pix);
}
}
