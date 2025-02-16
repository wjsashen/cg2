#include <fstream>
#include <vector>
#include <sstream>
#include <memory>
#include <cmath>
#include <limits>
#include <iostream>
#include <cmath>

//merge point and vector into one struct for easier calculation
struct Vec3d {
    float x, y, z;

    Vec3d(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}
    Vec3d operator*(double s) const {
        return Vec3d(x * s, y * s, z * s);
    }
    Vec3d operator/(double s) const {
        return Vec3d(x / s, y / s, z / s);
    }

    Vec3d operator+(const Vec3d& other) const {
        return Vec3d(x + other.x, y + other.y, z + other.z);
    }
    Vec3d operator-(const Vec3d& other) const {
        return Vec3d(x - other.x, y - other.y, z - other.z);
    }
    Vec3d cross(const Vec3d& b) const {
        return Vec3d(
            y * b.z - z * b.y,
            z * b.x - x * b.z,
            x * b.y - y * b.x
        );
    }
    //for debug
    // Overload << operator to print Vec3d objects
    friend std::ostream& operator<<(std::ostream& os, const Vec3d& vec) {
        os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
        return os;
    }
    float dot(const Vec3d& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    double lengthSquared() const {
        return x * x + y * y + z * z;
    }
    // Length (magnitude) of the vector
    float length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec3d norm() const {
        double len = sqrt(x * x + y * y + z * z);
        return (len > 0) ? Vec3d(x / len, y / len, z / len) : Vec3d(0, 0, 0);
    }
};

struct Ray {
    float x, y, z;
    float dx, dy, dz;
    Ray(const Vec3d& o, const Vec3d& d) 
        : x(o.x), y(o.y), z(o.z), dx(d.x), dy(d.y), dz(d.z) {}
};

struct Material {
    Vec3d diffuse, specular;
    float ka, kd, ks, shininess;
};

struct Light {
    Vec3d position;
    float intensity;
    bool isDirectional;
};
struct Camera {
    Vec3d viewDir, upDir;
    Vec3d eye;
    double vfov_degrees, w, h;
    double vfov_rad() const {
        return vfov_degrees * (M_PI / 180.0);
    }
};
double fHeight;
struct Color {
    double r, g, b;
    Color(double r = 0, double g = 0, double b = 0) : r(r), g(g), b(b) {}
    //for debug print
    friend std::ostream& operator<<(std::ostream& os, const Color& color) {
        os << "(" << color.r << ", " << color.g << ", " << color.b << ")";
        return os;
    }
};

//this is for expand the code to various types of shapes
struct Objects {
    virtual ~Objects() = default;
    virtual Color getColor() const = 0;
};
struct Sphere : public Objects {
    Vec3d center;
    double radius;
    Color color;
    Sphere(Vec3d c, double r, Color col) : center(c), radius(r), color(col) {}
    Color getColor() const override {
        return color;
    }
};
struct Cylinder : public Objects {
    Vec3d center;
    Vec3d dir;
    double radius,length;
    Color color;
    Cylinder(Vec3d c, Vec3d d,double r,double l, Color col) : center(c), dir(d),radius(r), length(l),color(col) {}
    Color getColor() const override {
        return color;
    }
};
struct Scene {
    Camera camera;
    Color bkgcolor, temp_mtlcolor;
    std::vector<std::unique_ptr<Objects>> objects;
};
// Intersection struct
struct Intersection {
    bool hit;
    Vec3d point, normal;
    Material material;
};

// Shade function using Blinn-Phong
Vec3d ShadeRay(const Intersection& inter, const Vec3d& viewDir, const std::vector<Light>& lights) {
    Vec3d color(0, 0, 0);
    for (const auto& light : lights) {
        Vec3d lightDir = light.isDirectional ? light.position.norm() : (light.position - inter.point).norm();
        Vec3d halfVec = (viewDir + lightDir).norm();
        float diff = std::max(inter.normal.dot(lightDir), 0.0f);
        float spec = std::pow(std::max(inter.normal.dot(halfVec), 0.0f), inter.material.shininess);
        Vec3d diffuseColor = inter.material.diffuse * inter.material.kd * diff;
        Vec3d specularColor = inter.material.specular * inter.material.ks * spec;
        color = color + (diffuseColor + specularColor) * light.intensity;
    }
    return color;
}

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

            } else if (type == "bkgcolor" || type == "mtlcolor") {
                double r, g, b;
                if (!(words >> r >> g >> b)) {
                    throw std::runtime_error("Expected 3 float values for " + type);
                }
                if (r < 0 || r > 1 || g < 0 || g > 1 || b < 0 || b > 1) {
                    throw std::runtime_error("Color values must be between 0 and 1");
                }
                Color color(r, g, b);
                if (type == "bkgcolor") scene.bkgcolor = color;
                else scene.temp_mtlcolor = color;

            } else if (type == "sphere") {
                float cx, cy, cz, radius;
                if (!(words >> cx >> cy >> cz >> radius)) {
                    throw std::runtime_error("Expected 4 float values for sphere");
                }
                //std::cout << "Sphere values: " << cx << " " << cy << " " << cz << " " << radius << std::endl;
                if (radius <= 0) {
                    throw std::runtime_error("Sphere radius must be positive");
                }
                scene.objects.push_back(std::unique_ptr<Sphere>(
                    new Sphere(Vec3d(cx, cy, cz), radius, Color(scene.temp_mtlcolor))));

            } else if (type == "parallel") {
                if (!(words >> fHeight)) {
                    throw std::runtime_error("Expected 1 float value for parallel");
                }
                if (fHeight <= 0) {
                    throw std::runtime_error("Height must be positive");
                }

            } else if (type == "cylinder") {
                float cx, cy, cz, dx, dy, dz, radius, length;
                if (!(words >> cx >> cy >> cz >> dx >> dy >> dz >> radius >> length)) {
                    throw std::runtime_error("Expected 8 float values for cylinder");
                }
                if (radius <= 0) {
                    throw std::runtime_error("Cylinder radius must be positive");
                }
                if (length <= 0) {
                    throw std::runtime_error("Cylinder length must be positive");
                }
                Vec3d dir(dx, dy, dz);
                if (dir.length() == 0) {
                    throw std::runtime_error("Cylinder direction vector cannot be zero");
                }
                dir = dir.norm();
                scene.objects.push_back(std::unique_ptr<Cylinder>(
                    new Cylinder(Vec3d(cx, cy, cz), dir, radius, length, Color(scene.temp_mtlcolor))));
            }

            // Check for extra unexpected parameters
            std::string extra;
            if (words >> extra) {
                throw std::runtime_error("Unexpected extra parameters");
            }

        } catch (const std::runtime_error& e) {
            std::cerr << "Error on line " << lineNumber << ": " << e.what() << std::endl;
            std::cerr << "Line content: " << line << std::endl;
            // You can choose to either continue parsing or exit here
            // return; // Uncomment to exit on first error
        }
    }

    // Validate that all required camera parameters were set
    if (scene.camera.w <= 0 || scene.camera.h <= 0) {
        throw std::runtime_error("Image size not set or invalid");
    }
    if (scene.camera.vfov_degrees <= 0) {
        throw std::runtime_error("Field of view not set or invalid");
    }
}
bool intersectRaySphere(const Ray& ray, const Sphere& sphere, float& t) {
    //quadratic approach for intersect
    float cx = sphere.center.x, cy = sphere.center.y, cz = sphere.center.z;
    float r = sphere.radius;

    //ray is normed so a is 1
    float a = ray.dx * ray.dx + ray.dy * ray.dy + ray.dz * ray.dz;
    float b = 2 * (ray.dx * (ray.x - cx) + ray.dy * (ray.y - cy) + ray.dz * (ray.dz - cz));
    float c = (ray.x - cx) * (ray.x - cx) + (ray.y - cy) * (ray.y - cy) + (ray.dz - cz) * (ray.dz - cz) - r * r;
    float d = b * b - 4 * a * c;
    if (d < 0) return false;
    float t1 = (-b - sqrt(d)) / (2 * a);
    float t2 = (-b + sqrt(d)) / (2 * a);
    t = (t1 > 0) ? t1 : (t2 > 0 ? t2 : -1);
    return t > 0;
}
bool intersectRayCylinder(const Ray &ray, const Cylinder &cylinder, float &t) {//for cylinder
  Vec3d rayo(ray.x, ray.y, ray.z);
  Vec3d rayDir(ray.dx, ray.dy, ray.dz);

    /*  
  perpendicular calculation  doesnt work well tho
  Vec3d co = rayo - cylinder.center;

  Vec3d d_perp = rayDir - cylinder.dir * rayDir.dot(cylinder.dir);
  Vec3d co_perp = co - cylinder.dir * co.dot(cylinder.dir);
  float A = d_perp.dot(d_perp);
  float B = 2 * d_perp.dot(co_perp);
  float C = co_perp.dot(co_perp) - cylinder.radius * cylinder.radius;*/

  //given cylinder equation  r² = x² + y² and cap is defined thru zmin<=z<=zmax
  // and go thru the deduction of putting ray equation to x,y,z
    // how does align with axis make things easier?
    float A = ray.dx * ray.dx + ray.dy * ray.dy;
    float B = 2*(ray.x*ray.dx+ray.y+ray.dy);
    float C = ray.x*ray.x+ray.y*ray.y-cylinder.radius*cylinder.radius;
    float D = B * B - 4 * A * C;
    if (D < 0 || fabs(A) < 1e-6) {
        //std::cout << "determinant is <0, no intersect" << std::endl;
        return false; 
    }
    float t1 = (-B - sqrt(D)) /(2*A);
    float t2 = (-B + sqrt(D)) /(2*A);
    float zmin = cylinder.center.z-cylinder.length/2.0;
    float zmax = cylinder.center.z+cylinder.length/2.0;

    //not sure here how to check cap
    t = (t1 > 0) ? t1 : (t2 > 0 ? t2 : -1);
    Vec3d z_hitv1 = rayo + rayDir *t1;
    float z_hit1= z_hitv1.z;
    Vec3d z_hitv2 = rayo + rayDir *t2;
    float z_hit2 = z_hitv2.z;
    bool res = (t>0) && ((z_hit1>=zmin&&z_hit1<=zmax)||(z_hit2>=zmin&&z_hit2<=zmax));
    return res;

}

Color paral(int x, int y, Vec3d orth, double fHeight, Scene& scene) { //for parallel project
    Ray r(scene.camera.eye, orth);
    Color pixelColor = scene.bkgcolor; // if it not hit anything
    double t_min = std::numeric_limits<double>::infinity();
    // Check for intersections with all objects in the scene
    for (const auto& obj : scene.objects) {
        Sphere* sphere = dynamic_cast<Sphere*>(obj.get());
        if (sphere) {
            float t=fHeight;//dist along the ray
            if (intersectRaySphere(r, *sphere, t) && t < t_min) {
                t_min = t;  //make sure nearest intersect get the color
                pixelColor = sphere->getColor();  
                //std::cout << "Intersected with Sphere at t = " << t << ", Color: " << pixelColor << std::endl;
            }
        }
    }
    return pixelColor;
}

Color trace(int x,int y,Vec3d ul, Vec3d delta_h, Vec3d delta_v, Scene& scene){  //for pespective projection ray trace
    /* this is the fov based approach to calculate the ray dir, which gave a slightly different result in shapes
    double fov = tan(scene.camera.vfov_rad()/2);
    double vx = (2*(x+0.5)/width-1)*ar*fov;
    double vy =(1-2*(y+0.5)/height)*fov;
    Vec3d viewDir = scene.camera.viewDir.norm();
    Vec3d right = viewDir.cross(scene.camera.upDir).norm();  
    Vec3d rayDir = (viewDir + right * vx + scene.camera.upDir * vy).norm();  
    */

    // (corner-based):
    // ray target position on the view window
    Vec3d vwPosition = ul + delta_h * x + delta_v * y;
    Vec3d rayDir = (vwPosition - scene.camera.eye).norm();
    Ray r = Ray(scene.camera.eye,rayDir);
    float t_min = std::numeric_limits<float>::max();
    Color pixelColor = scene.bkgcolor; // if it not hit anything
    // when representing in calculation
    // ray = rayo + scale*rayDir. 
    // scale=0 we at rayo when 1 we at the window(focal) when above beyond window
    for (const auto &obj : scene.objects) {
      // Check for intersection with Spheres one by one
      Sphere *sphere = dynamic_cast<Sphere *>(obj.get());
      if (sphere) {
        float t;
        if (intersectRaySphere(r, *sphere, t) && t < t_min) {
          t_min = t;
          pixelColor = sphere->getColor();
        }
      }
      // Check for intersection with Cylinder, check logic here.
      Cylinder *cylinder = dynamic_cast<Cylinder *>(obj.get());
      if (cylinder) {
        float t;
        if (intersectRayCylinder(r, *cylinder, t) && t < t_min) {
          t_min = t;
          pixelColor = cylinder->getColor();
        }
      }
    }

    return pixelColor;
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
    std::vector<std::vector<Color>> image2(height, std::vector<Color>(width, Color(0, 0, 0)));
    
    std::string perspective_filename = basename + "_perspective.ppm";
    std::string parallel_filename = basename + "_parallel.ppm";
    
    std::ofstream output(perspective_filename);
    output << "P3\n" << width << " " << height << "\n255\n";  
    std::ofstream output2(parallel_filename);
    output2 << "P3\n" << width << " " << height << "\n255\n";  

    //the view window init is here, (maybe need refactor it for later assignments)
    //didn't seperate it to a func due to easier param passing
    // and it will be modified by the projection ways just easier for checking:
    
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


    //double fWidth = fHeight * ar;
    // Compute the four corners of the view window for parallel projection by chaning the param of 4 corners(only need 3 corners in code)
    //ray not converge, so directly use eye as "center" in corners cal(d=0)
    ul = scene.camera.eye - right * (vw * 0.5) + up * (vh * 0.5);;//focal is 0
    ur = scene.camera.eye + right * (vw * 0.5) + up * (vh * 0.5); ;
    ll = scene.camera.eye + right * (vw * 0.5) + up * (vh * 0.5); ;
    //Vec3d lr = scene.camera.eye - up * (fHeight / 2.0) + right * (fWidth / 2.0);
    Vec3d orthogonalDir = right.cross(up).norm(); //since ray dir will be same for paral
    // Trace for parallel projection image (second render)
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            image2[y][x] = paral(x, y, orthogonalDir, fHeight, scene);
            output2 << (int)(image2[y][x].r * 255) << " "
                    << (int)(image2[y][x].g * 255) << " "
                    << (int)(image2[y][x].b * 255) << " ";
        }
        output2 << "\n";
    }
    output2.close(); 
    std::cout << "Rendering complete. Images saved as 'output_perspective.ppm' and 'output_parallel.ppm'." << std::endl;
    return 0;
}
