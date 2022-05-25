#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "vector.cpp"
#include <vector>
#include <math.h>
#include <iostream>
#include <limits>
#include <cmath>
#include <chrono>

int main(){
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

    Vector red = Vector(0.9,0,0);
    Vector green = Vector(0,0.9,0);
    Vector blue = Vector(0,0,0.9);
    Vector black = Vector(0,0,0);
    
    std::vector<Sphere> objects = {Sphere(Vector(0,1000,0),940,red, false), 
                                   Sphere(Vector(0,0,-1000),940,green, false),
                                   Sphere(Vector(0,-1000,0),990,blue, false),
                                   Sphere(Vector(0,0,1000),940, blue + red, false),
                                   Sphere(Vector(1000, 0, 0), 940, blue + green, false),
                                   Sphere(Vector(-1000, 0, 0), 940, red + green, false),
                                   Sphere(Vector(0,0,0),10, blue + red + green, false)};
    
    //std::vector<Sphere> objects = {Sphere(Vector(0,0,0),10, blue + red + green)};
    Scene scene = Scene(objects);
    Vector camera = Vector(0,0,55);
    PointLight light_source = PointLight(Vector(-10,20,40), 2E10);

    int W = 512;
	int H = 512;
    double alpha = M_PI/3;

    int ray_depth = 5;
    int number_of_paths = 1000;
    //int antialiasing_sample_size = 10;

	std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for schedule(dynamic , 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
            Vector direction = Vector(j + 0.5 - W/2, (H - i - 1) + 0.5 - H/2, -W/(2*tan(alpha/2)));
            direction.normalize();
            Ray ray = Ray(camera,  direction);
            Intersection closest = Intersection(false,0,black);
            int n_objs = scene.objects.size();
            for (int k = 0; k < n_objs; k++) {
                scene.objects[k].intersect(&ray, &closest);
            }

            if (closest.does_int){
                Vector P = camera + direction*closest.int_time;
                Vector normal_vector = P - closest.sphere_center;
                normal_vector.normalize();
                double epsilon = 0.0000001;

                Vector pixel_color = Vector(0,0,0);
                for (int i = 0; i < number_of_paths; i += 1){
                    /*
                    Vector rand = Vector(0,0,0);
                    boxMuller(1, rand[0],rand[1]);
                    Vector random_dir = direction + rand;
                    Ray rand_ray(camera,random_dir);

                    Intersection rand_intersection = Intersection(false,0,black);
                    int n_objs = scene.objects.size();
                    for (int k = 0; k < n_objs; k++) {
                        scene.objects[k].intersect(&rand_ray, &rand_intersection);
                    }
                    if (rand_intersection.does_int) {
                        Vector P_rand = camera + random_dir*rand_intersection.int_time;
                        Vector normal_vector_rand = P_rand - rand_intersection.sphere_center;
                        normal_vector_rand.normalize();
                        double epsilon = 0.0000001;
                        pixel_color = pixel_color + light_source.getColor(P_rand + normal_vector_rand*epsilon , scene, rand_intersection.color, normal_vector_rand, rand_intersection.is_mirror, rand_ray, ray_depth);
                    }
                    */
                    pixel_color = pixel_color + light_source.getColor(P + normal_vector*epsilon , scene, closest.color, normal_vector, closest.is_mirror, ray, ray_depth);
                }
                pixel_color = pixel_color/(number_of_paths);

                image[(i * W + j) * 3 + 0] = std::min(255, (int) pow(pixel_color[0], 1./2.2));
                image[(i * W + j) * 3 + 1] = std::min(255, (int) pow(pixel_color[1], 1./2.2));
                image[(i * W + j) * 3 + 2] = std::min(255, (int) pow(pixel_color[2], 1./2.2));
            }
		}
	}
	stbi_write_png("spheres.png", W, H, 3, &image[0], 0);

    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cout << "Rendering time:" << ms_int.count() << "ms" << std::endl;



	return 0;
}