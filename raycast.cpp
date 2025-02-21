#include <iostream>
#include "raycast.h"
#include "parse.h"
#include "raycast.h"
#include <fstream>
#include <vector>
#include <sstream>
#include <memory>
#include <cmath>
#include <limits>
#include <iostream>
#include <cmath>


    
bool intersectRaySphere(const Ray& ray, const Sphere& sphere, float& t) {
    //intersect func
        float cx = sphere.center.x, cy = sphere.center.y, cz = sphere.center.z;
        float r = sphere.radius;
        //fixed a issue from 1a
        float a = ray.dx * ray.dx + ray.dy * ray.dy + ray.dz * ray.dz;
        float b = 2 * (ray.dx * (ray.x - cx) + ray.dy * (ray.y - cy) + ray.dz * (ray.z - cz));
        float c = (ray.x - cx) * (ray.x - cx) + (ray.y - cy) * (ray.y - cy) + (ray.z - cz) * (ray.z - cz) - r * r;
    
        float d = b * b - 4 * a * c;
        if (d < 0) return false; 
    
        float t1 = (-b - sqrt(d)) / (2 * a);
        float t2 = (-b + sqrt(d)) / (2 * a);
    
        if (t1 > 0) {
            t = t1;
        } else if (t2 > 0) {
            t = t2;
        } else {
            return false; 
        }
        return true;
    }
    

bool isInShadow(const Vec3d& point, const Vec3d& L, const Scene& scene) {
    //check if the point is in shadow area
        Ray shadowRay(point + L * 0.001, L); // offset to avoid self-intersection
    
        for (const auto& obj : scene.objects) {
            Sphere* sphere = dynamic_cast<Sphere*>(obj.get());
            if (sphere) {
                float t;
                if (intersectRaySphere(shadowRay, *sphere, t) && t > 0.001) {
                    return true; // Point is in shadow
                }
            }
        }
        return false; //light reaches point
    }
    

    Color shade(Scene& scene, float t_min, Vec3d rayDir,float t, Sphere sphere){
        //drag shade out to sperate func
        Color finalColor = scene.bkgcolor;
        Vec3d point,N;
        MaterialColor mt;
    
        point = scene.camera.eye + rayDir * t;
        mt = sphere.getColor();
    
        N = (point - sphere.center).norm();
    
        if (t_min < std::numeric_limits<float>::max()) {
            Vec3d viewDir = (rayDir * -1).norm();
            Color ambient = mt.color * mt.ka; 
        
            Color totalDiffuse(0, 0, 0);
            Color totalSpecular(0, 0, 0);
        
            // multiple lights in one ray
            for (const auto& light : scene.lights) {
                //TODO:distinguish light is dir/point
                Vec3d L = light->isPoint ?   (light->positionOrdir - point).norm():light->positionOrdir.norm()*-1;
                float dis = (light->positionOrdir - point).length(); // distance to light
                //with the default set, fatt is 1 when c1 c2 c3 is not in light params 
                float fatt = 1.0f / (light->c1 + light->c2 * dis + light->c3 * dis * dis);
                Vec3d H = (L + viewDir).norm();
                float df = std::max(N.dot(L), 0.0f);
                
                Color diffuse = mt.color * mt.kd * df;
                float sf = std::pow(std::max(N.dot(H), 0.0f), mt.shininess);
                Color specular = mt.specular * mt.ks * sf;
        
                float Si = isInShadow(point, L, scene) ? 0.0f : 1.0f;
    
                totalDiffuse = totalDiffuse+ diffuse * light->intensity * Si*fatt;
                totalSpecular = totalSpecular+specular * light->intensity * Si*fatt;
            }
        
            // Final color calculation
            finalColor = ambient + totalDiffuse + totalSpecular;
        }
        //clamp made the color to 0-1 value, but not seems neccessary some test cases with given color value
        //std::cout<<finalColor<<"final color before clamp is"<<std::endl;
        return finalColor.clamp();
    }

Color trace(int x, int y, Vec3d ul, Vec3d delta_h, Vec3d delta_v, Scene& scene) {
    
        Vec3d vwPosition = ul + delta_h * x + delta_v * y;
        Vec3d rayDir = (vwPosition - scene.camera.eye).norm();
        Ray r(scene.camera.eye, rayDir);
        float t_min = std::numeric_limits<float>::max();
        Color res = scene.bkgcolor;
        for (const auto &obj : scene.objects) {
          Sphere *sphere = dynamic_cast<Sphere *>(obj.get());
          if (sphere) {
            float t;
            if (intersectRaySphere(r, *sphere, t) && t < t_min) {
              t_min = t;
              //std::cout << "Hit sphere at t = " << t << " with center " << sphere->center << std::endl;
              res = shade(scene,t_min,rayDir,t,*sphere);
            }
          }
        }
        return res;
    }

bool isValidTxtFile(const std::string& filename) { //check file format
    return filename.size() >= 4 && filename.substr(filename.size() - 4) == ".txt";
} 

int main() {
    Scene scene;
    std::string filename;
    std::cout << "Enter the input file name: ";
    std::cin >> filename;
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Error: Could not open file '" << filename << "'." << std::endl;
        return 1; 
    }else if (!isValidTxtFile(filename)) {
        std::cerr << "Error: Invalid file format. Please provide a .txt file." << std::endl;
        return 1;
    }
    std::string basename = filename;
    size_t lastdot = basename.find_last_of(".");
    if (lastdot != std::string::npos) {
        basename = basename.substr(0, lastdot);
    }

    parse(filename, scene);
    int width = (int)scene.camera.w;
    int height = (int)scene.camera.h;
    double ar = (double)width / height;
    std::vector<std::vector<Color>> image(height, std::vector<Color>(width, Color(0, 0, 0)));
    std::string perspective_filename = basename + "_perspective.ppm";
    std::ofstream output(perspective_filename);
    output << "P3\n" << width << " " << height << "\n255\n";  


    //the view window init is here,
    // TODO: need refactor it for a sperate func but need a new struct to store return value, like a vector of Vec3d
    
    double vh = 2.0 * tan(scene.camera.vfov_rad() / 2);  //view window width and height in 3d world coord
    double vw = vh * ar; //the given w&h is in pixel so it need to change to view window measurement
    Vec3d viewDir = scene.camera.viewDir.norm();
    //pay attention to the cross sequence viewD first(right hand)
    Vec3d right = viewDir.cross(scene.camera.upDir).norm();//(u in lecture notes, the orthogonal to the plane)
    // (up is the v in lecture note, the norm of orth of u&viewD )
    Vec3d up = right.cross(viewDir).norm(); //v in lecture notes


    // four corners of view window
    //d*n in notes
    Vec3d center = scene.camera.eye + viewDir; // d(focal distance) is set to arbitary 1, that's where frustum narrow to the center of vw
    Vec3d ul = center - right * (vw * 0.5) + up * (vh * 0.5); // Upper-left corner
    Vec3d ur = center + right * (vw * 0.5) + up * (vh * 0.5); // Upper-right corner
    Vec3d ll = center - right * (vw * 0.5) - up * (vh * 0.5); // Lower-left corner
    // step vectors
    Vec3d delta_h = (ur - ul) / width;
    Vec3d delta_v = (ll - ul) / height;
    // Trace for perspective image (first render)
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            image[y][x] = trace(x, y, ul, delta_h, delta_v, scene);
            //or no need for buffer here whatever
            output << (int)(image[y][x].r * 255) << " "
                   << (int)(image[y][x].g * 255) << " "
                   << (int)(image[y][x].b * 255) << " ";
        }
        output << "\n";
    }
    output.close();  


  
    std::cout << "Rendering complete. Images saved as 'output_perspective.ppm' and 'output_parallel.ppm'." << std::endl;
    return 0;
}
