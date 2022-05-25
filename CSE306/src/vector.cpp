#include <vector>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <random>

static std::default_random_engine engine (10);
static std::uniform_real_distribution<double> uniform(0,1);

double pi = M_PI;

class Vector {
    private:
        double coords[3];
    public:
        explicit Vector(double x = 0., double y = 0., double z = 0.){
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        };
        double norm2() const {
		    return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
	    }
	    double norm() const {
		    return sqrt(norm2());
	    }
	    void normalize() {
		    double n = norm();
		    coords[0] /= n;
		    coords[1] /= n;
		    coords[2] /= n;
	    }

        double operator[](int i) const { return coords[i]; };
	    double& operator[](int i) { return coords[i]; };
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}


class Ray {
    public:
        Vector origin;
        Vector direction;
        Ray(Vector O, Vector u) {
            origin = O;
            direction = u;
        }
};




class Intersection {
    public:
        bool does_int;
        double int_time;
        Vector color;
        Vector sphere_center;
        bool is_mirror;

        Intersection(bool b, double t, Vector clr){
            color = clr;
            does_int = b;
            if (b == true) {
                int_time = t;
            }
            else {
                int_time = -1;
            }
        }
};

class Sphere {
    public:
        Vector center;
        double radius;
        Vector albedo;
        bool is_mirror;

        Sphere(){
            center = Vector(0,0,0);
            radius = 0;
            albedo = Vector(0,0,0);
        }
        Sphere(Vector c, double r, Vector color, bool mirror){
            center = c;
            radius = r;
            albedo = color;
            is_mirror = mirror;
        }
        bool intersect(Ray* ray, Intersection* int_data) {
            Vector O = ray->origin;
            Vector u = ray->direction;
            Vector n = O - center;
            double delta = (pow(dot(u,O - center),2) - (n.norm2()-pow(radius,2)));
            //std::cout << delta << std::endl;
            if (delta < 0) {
                return false;
            }else {
                double t1 = dot(u, center - O) - sqrt(delta);
                double t2 = dot(u, center - O) + sqrt(delta);
                if (t2 < 0){
                    return false;
                }
                else if (t1 < 0) {
                    if (!(int_data->does_int)) {
                        int_data->does_int = true;
                        int_data->int_time = t2;
                        int_data->color = albedo;
                        int_data->sphere_center = center;
                        int_data->is_mirror = is_mirror;
                        return true;
                    }
                    else {
                        if (t2 < int_data->int_time) {
                            int_data->int_time = t2;
                            int_data->color = albedo;
                            int_data->sphere_center = center;
                            int_data->is_mirror = is_mirror;
                            return true;
                        }
                        else {
                            return false;
                        }
                    }
                }
                else{
                    if (!(int_data->does_int)) {
                        int_data->does_int = true;
                        int_data->int_time = t1;
                        int_data->color = albedo;
                        int_data->sphere_center = center;
                        int_data->is_mirror = is_mirror;
                        return true;
                    }
                    else {
                        if (t1 < int_data->int_time) {
                            int_data->int_time = t1;
                            int_data->color = albedo;
                            int_data->sphere_center = center;
                            int_data->is_mirror = is_mirror;
                            return true;
                        }
                        else {
                            return false;
                        }
                    }
                    
                }
            }
        }
};

class Scene {
    public:
        std::vector<Sphere> objects;
    
        Scene(std::vector<Sphere> object_vect) {
            objects = object_vect;
        }
};

void boxMuller(double stdev, double &x, double &y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2*log(r1))*cos(2*pi*r2)*stdev;
    y = sqrt(-2*log(r1))*sin(2*pi*r2)*stdev;
} 


class PointLight {
    public:
        Vector origin;
        double intensity;
        PointLight(Vector O, double I){
            origin = O;
            intensity = I;
        }

        Vector random_cos(Vector N){
            double r1 = uniform(engine);
            double r2 = uniform(engine);

            double x = cos(2*pi*r1)*sqrt(1-r2);
            double y = sin(2*pi*r1)*sqrt(1-r2);
            double z = sqrt(r2);

            Vector T1 = Vector(-N[1], N[0], 0);
            T1.normalize();
            Vector T2 = cross(N, T1);
            return x*T1 + y*T2 + z*N;
        }

        //returns color values
        Vector getColor(Vector P, Scene scene, Vector color, Vector normal, bool is_mirror, Ray source_ray, int ray_depth){
            if (ray_depth < 0) {
                return Vector(0,0,0);
            }

            if (is_mirror) {
                Vector reflection_direction = source_ray.direction - 2*dot(source_ray.direction, normal)*normal;
                //reflection_direction.normalize();
                Ray reflection = Ray(P, reflection_direction);

                Intersection closest = Intersection(false,0,Vector(0,0,0));
                int n_objs = scene.objects.size();
                for (int k = 0; k < n_objs; k++) {
                    scene.objects[k].intersect(&reflection, &closest);
                }

                if (closest.does_int){
                    Vector P_next = P + reflection_direction*closest.int_time;
                    Vector normal_next = P_next - closest.sphere_center;
                    normal_next.normalize();
                    double epsilon = 0.0000001;
                    return this->getColor(P_next + normal_next*epsilon, scene, closest.color, normal_next ,closest.is_mirror, reflection, ray_depth);
                }
                else {
                    return Vector(0,0,0);
                }
            }

            Vector direction = origin - P;
            direction.normalize();
            Vector dist = origin - P;
            double distance = dist.norm();
            Ray to_light_source = Ray(P, direction);

            Intersection closest = Intersection(false,0,Vector(0,0,0));
            int n_objs = scene.objects.size();
            for (int k = 0; k < n_objs; k++) {
                scene.objects[k].intersect(&to_light_source, &closest);
            }

            Vector Lo = Vector(0,0,0);
            if (closest.does_int){
                if (closest.int_time >= distance){
                    Lo = Lo + (intensity/(4*pi*pow(distance,2)))*(color/pi)*(std::max(dot(normal,direction),0.));
                }
                Ray random_ray = Ray(P, random_cos(normal));

                Intersection next_intersection = Intersection(false,0,Vector(0,0,0));
                int n_objs = scene.objects.size();
                for (int k = 0; k < n_objs; k++) {
                    scene.objects[k].intersect(&random_ray, &next_intersection);
                }
                if (next_intersection.does_int) {
                    Vector P_next = P + next_intersection.int_time*random_ray.direction;
                    Vector normal_next = P_next - next_intersection.sphere_center;
                    normal_next.normalize();
                    double epsilon = 0.0000001;
                    Vector computed = getColor(P_next + normal_next*epsilon, scene, next_intersection.color, normal_next, next_intersection.is_mirror, random_ray, ray_depth-1);
                    computed[0] *= color[0];
                    computed[1] *= color[1];
                    computed[2] *= color[2];
                    Lo = Lo + computed;
                    return Lo;
                }
                else {
                    return Lo;
                }
                 
            }
            else {
                Lo = Lo + (intensity/(4*pi*pow(distance,2)))*(color/pi)*(std::max(dot(normal,direction),0.));
                
                Ray random_ray = Ray(P, random_cos(normal));

                Intersection next_intersection = Intersection(false,0,Vector(0,0,0));
                int n_objs = scene.objects.size();
                for (int k = 0; k < n_objs; k++) {
                    scene.objects[k].intersect(&random_ray, &next_intersection);
                }
                if (next_intersection.does_int) {
                    Vector P_next = P + next_intersection.int_time*random_ray.direction;
                    Vector normal_next = P_next - next_intersection.sphere_center;
                    normal_next.normalize();
                    double epsilon = 0.0000001;
                    Vector computed = getColor(P_next + normal_next*epsilon, scene, next_intersection.color, normal_next, next_intersection.is_mirror, random_ray, ray_depth-1);
                    computed[0] *= color[0];
                    computed[1] *= color[1];
                    computed[2] *= color[2];
                    Lo = Lo + computed;
                    return Lo;
                }
                else {
                    return Lo;
                }
            }
        }
};