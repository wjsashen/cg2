#include "parse.h"

#include <fstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <memory>
#include <cmath>
#include <limits>
#include <iostream>
#include <cmath>


bool isWhitespaceOnly(const std::string& str) {//for parse check
    return str.find_first_not_of(" \t\r\n") == std::string::npos;
}
void parse(const std::string& filename, Scene& scene) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::string type;
    int lineNumber = 0;

    while (std::getline(file, line)) {
        lineNumber++;
        if (isWhitespaceOnly(line)) {
            continue;
        }
        std::istringstream words(line);
        words >> type;

        try {
            if (type == "eye" || type == "viewdir" || type == "updir") {
                float x, y, z;
                if (!(words >> x >> y >> z)) {
                    throw std::runtime_error("Expected 3 float values for " + type);
                }
                Vec3d vec(x, y, z);
                if (type == "eye") scene.camera.eye = vec;
                else if (type == "viewdir") scene.camera.viewDir = vec;
                else scene.camera.upDir = vec;

            } else if (type == "vfov") {
                if (!(words >> scene.camera.vfov_degrees)) {
                    throw std::runtime_error("Expected 1 float value for vfov");
                }
                if (scene.camera.vfov_degrees <= 0 || scene.camera.vfov_degrees >= 180) {
                    throw std::runtime_error("vfov must be between 0 and 180 degrees");
                }

            } else if (type == "imsize") {
                if (!(words >> scene.camera.w >> scene.camera.h)) {
                    throw std::runtime_error("Expected 2 integer values for imsize");
                }
                if (scene.camera.w <= 0 || scene.camera.h <= 0) {
                    throw std::runtime_error("Image dimensions must be positive");
                }

            } else if (type == "bkgcolor") {
                double r, g, b;
                
                if (!(words >> r >> g >> b)) {
                    throw std::runtime_error("Expected 3 float values for bkgcolor");
                }
                if (r < 0 || r > 1 || g < 0 || g > 1 || b < 0 || b > 1) {
                    throw std::runtime_error("Color values must be between 0 and 1");
                }
                
                scene.bkgcolor = Color(r, g, b);
            } 
            else if (type == "mtlcolor") {
                double dr, dg, db, sr, sg, sb, ka, kd, ks, shininess;
            
                if (!(words >> dr >> dg >> db >> sr >> sg >> sb >> ka >> kd >> ks >> shininess)) {
                    std::cerr << "Failed to parse mtlcolor. Extracted values: " 
                              << dr << " " << dg << " " << db << " " << sr << " " << sg << " " 
                              << sb << " " << ka << " " << kd << " " << ks << " " << shininess 
                              << std::endl;
                    throw std::runtime_error("Expected 10 values for mtlcolor (dr dg db sr sg sb ka kd ks shininess)");
                }
                
                // Validate diffuse and specular color values
                if ((dr < 0 || dr > 1) || (dg < 0 || dg > 1) || (db < 0 || db > 1) ||
                    (sr < 0 || sr > 1) || (sg < 0 || sg > 1) || (sb < 0 || sb > 1)) {
                    throw std::runtime_error("Diffuse and specular color values must be between 0 and 1");
                }
            
                // Validate coefficients and shininess
                if (ka < 0 || kd < 0 || ks < 0 || shininess < 0) {
                    throw std::runtime_error("Material coefficients (ka, kd, ks) and shininess must be non-negative");
                }
            
                // Construct MaterialColor and assign to scene
                scene.temp_mtlcolor = MaterialColor(
                    Color(dr, dg, db),    // Diffuse color
                    Color(sr, sg, sb),    // Specular color
                    ka, kd, ks, shininess // Material properties
                );
            }
            else if (type == "sphere") {
                float cx, cy, cz, radius;
                if (!(words >> cx >> cy >> cz >> radius)) {
                    throw std::runtime_error("Expected 4 float values for sphere");
                }
                //std::cout << "Sphere values: " << cx << " " << cy << " " << cz << " " << radius << std::endl;
                if (radius <= 0) {
                    throw std::runtime_error("Sphere radius must be positive");
                }
                scene.objects.push_back(std::unique_ptr<Sphere>(
                    new Sphere(Vec3d(cx, cy, cz), radius, MaterialColor(scene.temp_mtlcolor))));

            } else if (type == "light") {
                float x, y, z;
                int isPoint;
                float intensity;
            
                if (!(words >> x >> y >> z >> isPoint >> intensity)) {
                    throw std::runtime_error("Expected format: light x y z isPoint intensity");
                }
                
                if (intensity < 0) {
                    throw std::runtime_error("Light intensity must be non-negative");
                }
            
                // Debugging
                //std::cout << "Parsed light: Position(" << x << ", " << y << ", " << z 
                  //        << ") isPoint: " << isPoint << " Intensity: " << intensity << std::endl;
                
                float c1 = 1.0f, c2 = 0.0f, c3 = 0.0f; //  no attenuation
                if (words >> c1 >> c2 >> c3) {
                    std::cout<<"have att params"<<std::endl;
                }
                Light light;
                light.c1 = c1;
                light.c2 =c2;
                light.c3=c3;
                light.positionOrdir = Vec3d(x, y, z);
                light.isPoint = (isPoint != 0); // Convert int to boolean
                light.intensity = intensity;
                scene.lights.push_back(std::make_unique<Light>(light)); 
            }
            
            std::string extra;
            if (words >> extra) {
                throw std::runtime_error("Unexpected extra parameters");
            }

        } catch (const std::runtime_error& e) {
            std::cerr << "Error on line " << lineNumber << ": " << e.what() << std::endl;
            std::cerr << "Line content: " << line << std::endl;
            // return; // Uncomment to exit on first error
        }
    }

    if (scene.camera.w <= 0 || scene.camera.h <= 0) {
        throw std::runtime_error("Image size not set or invalid");
    }
    if (scene.camera.vfov_degrees <= 0) {
        throw std::runtime_error("Field of view not set or invalid");
    }
}