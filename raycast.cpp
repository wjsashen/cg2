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
    Vec3d reflect(const Vec3d& normal) const {
        return *this - normal * (2 * this->dot(normal));
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
    // Direct component access - these are now public member variables, not methods
    float x, y, z;        // origin components
    float dx, dy, dz;     // direction components
    
    // Constructor remains the same
    Ray(const Vec3d& o, const Vec3d& d)
        : x(o.x), y(o.y), z(o.z), 
          dx(d.x), dy(d.y), dz(d.z) {}
          
    // Optional: Add method to get vectors if needed
    Vec3d getOrigin() const { return Vec3d(x, y, z); }
    Vec3d getDirection() const { return Vec3d(dx, dy, dz); }
};


struct Light {
    Vec3d position;
    bool isPoint;
    float intensity;
};
struct Camera {
    Vec3d viewDir, upDir;
    Vec3d eye;
    double vfov_degrees, w, h;
    double vfov_rad() const {
        return vfov_degrees * (M_PI / 180.0);
    }
};
struct Color {
    double r, g, b;
    Color(double r = 0, double g = 0, double b = 0) : r(r), g(g), b(b) {}
    
    Color operator*(double s) const { return Color(r * s, g * s, b * s); }
    Color operator+(const Color& other) const { return Color(r + other.r, g + other.g, b + other.b); }

    // Clamp color values to [0,1]
    Color clamp() const {
        return Color(
            std::max(0.0, std::min(1.0, r)),
            std::max(0.0, std::min(1.0, g)),
            std::max(0.0, std::min(1.0, b))
        );
    }
    friend std::ostream& operator<<(std::ostream& os, const Color& color) {
        os << "(" << color.r << ", " << color.g << ", " << color.b << ")";
        return os;
    }
};

struct MaterialColor {
    Color color;        // Base color
    Color specular;
    double ka, kd, ks;   
    double shininess;

    MaterialColor(
        const Color& c = Color(), 
        const Color& s = Color(), 
        double ka = 0, double kd = 0, double ks = 0, double shininess = 0
    ) : color(c),  specular(s), ka(ka), kd(kd), ks(ks), shininess(shininess) {}

};
//this is for expand the code to various types of shapes
struct Objects {
    virtual ~Objects() = default;
    virtual MaterialColor getColor() const = 0;
    Vec3d center;
    Vec3d getCenter(){
        return center;
    }
};
struct Sphere : public Objects {
    double radius;
    MaterialColor material; 
    Vec3d center;

    Sphere(Vec3d c, double r, MaterialColor mat)  // Accept MaterialColor
        : center(c), radius(r), material(mat) {}

    MaterialColor getColor() const override {
        return material;  // Return MaterialColor
    }

};

struct Cylinder : public Objects {
    Vec3d dir;
    Vec3d center;

    double radius,length;
    MaterialColor color;
    Cylinder(Vec3d c, Vec3d d,double r,double l, Color col) : center(c), dir(d),radius(r), length(l),color(col) {}
    MaterialColor getColor() const override {
        return color;
    }
};
struct Scene {
    Camera camera;
    Color bkgcolor;
    std::vector<std::unique_ptr<Light>> lights;
    MaterialColor temp_mtlcolor;
    std::vector<std::unique_ptr<Objects>> objects;
};
// Intersection struct
struct Intersection {
    bool hit;
    Vec3d point, normal;
    MaterialColor material;
};

class Plane : public Objects {
    public:
        Vec3d point;    // A point on the plane
        Vec3d normal;   // Plane normal vector
        MaterialColor material;
    
        Plane(const Vec3d& p, const Vec3d& n, const MaterialColor& mat) 
            : point(p), normal(n.norm()), material(mat) {}
    
        MaterialColor getColor() const override {
            return material;
        }
    };

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
                std::cout << "Parsed light: Position(" << x << ", " << y << ", " << z 
                          << ") isPoint: " << isPoint << " Intensity: " << intensity << std::endl;
            
                Light light;
                light.position = Vec3d(x, y, z);
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

    
    /*bool intersectRayPlane(const Ray& ray, const Plane& plane, float& t) {
        float denom = plane.normal.dot(ray.getDirection());
        
        if (std::abs(denom) > 1e-6) {  // Avoid division by zero
            Vec3d p0l0 = plane.point - ray.getOrigin();
            t = p0l0.dot(plane.normal) / denom;
            return (t >= 1e-6);  // Ensure intersection happens in front of the ray
        }
        return false;
    }*/
    
    bool intersectRaySphere(const Ray& ray, const Sphere& sphere, float& t) {
        float cx = sphere.center.x, cy = sphere.center.y, cz = sphere.center.z;
        float r = sphere.radius;
    
        // Corrected calculations
        float a = ray.dx * ray.dx + ray.dy * ray.dy + ray.dz * ray.dz;
        float b = 2 * (ray.dx * (ray.x - cx) + ray.dy * (ray.y - cy) + ray.dz * (ray.z - cz));
        float c = (ray.x - cx) * (ray.x - cx) + (ray.y - cy) * (ray.y - cy) + (ray.z - cz) * (ray.z - cz) - r * r;
    
        float d = b * b - 4 * a * c;
        if (d < 0) return false; // No real roots, no intersection
    
        // Compute t values
        float t1 = (-b - sqrt(d)) / (2 * a);
        float t2 = (-b + sqrt(d)) / (2 * a);
    
        // Fix: Ensure the correct `t` is chosen
        if (t1 > 0) {
            t = t1;
        } else if (t2 > 0) {
            t = t2;
        } else {
            return false; // Both intersections are behind the ray
        }
    
        return true;
    }
    

    bool isInShadow(const Vec3d& point, const Light& light, const Scene& scene) {
        Vec3d lightDir = (light.position - point).norm();
        Ray shadowRay(point + lightDir * 0.001, lightDir); // Offset to avoid self-intersection
    
        for (const auto& obj : scene.objects) {
            Sphere* sphere = dynamic_cast<Sphere*>(obj.get());
            if (sphere) {
                float t;
                if (intersectRaySphere(shadowRay, *sphere, t) && t > 0.001) {
                    return true; // Point is in shadow
                }
            }
        }
        return false; // No obstruction, light reaches point
    }
    



    Color trace(int x, int y, Vec3d ul, Vec3d delta_h, Vec3d delta_v, Scene& scene) {
    
        Vec3d vwPosition = ul + delta_h * x + delta_v * y;
        Vec3d rayDir = (vwPosition - scene.camera.eye).norm();
        Ray r(scene.camera.eye, rayDir);
        float t_min = std::numeric_limits<float>::max();
        Color finalColor = scene.bkgcolor;
        MaterialColor mt;
        Vec3d point, N;
        for (const auto &obj : scene.objects) {
          Sphere *sphere = dynamic_cast<Sphere *>(obj.get());
          if (sphere) {
            float t;
            if (intersectRaySphere(r, *sphere, t) && t < t_min) {
              t_min = t;
              point = scene.camera.eye + rayDir * t;
              mt = sphere->getColor();
              N = (point - sphere->center).norm();
              //std::cout << "Hit sphere at t = " << t << " with center " << sphere->center << std::endl;
                          }
          }
        }
        if (t_min < std::numeric_limits<float>::max()) {
            Vec3d viewDir = (rayDir * -1).norm();
            Color ambient = mt.color * mt.ka; 
        
            Color totalDiffuse(0, 0, 0);
            Color totalSpecular(0, 0, 0);
        
            // Iterate over ALL lights in the scene
            for (const auto& light : scene.lights) {
                Vec3d L = (light->position - point).norm();
                Vec3d H = (L + viewDir).norm();
                float df = std::max(N.dot(L), 0.0f);
                
                Color diffuse = mt.color * mt.kd * df;
                float sf = std::pow(std::max(N.dot(H), 0.0f), mt.shininess);
                Color specular = mt.specular * mt.ks * sf;
        
                // Check if the point is in shadow for this light
                float Si = isInShadow(point, *light, scene) ? 0.0f : 1.0f;
        
                // Add this light's contribution
                totalDiffuse = totalDiffuse+ diffuse * light->intensity * Si;
                totalSpecular = totalSpecular+specular * light->intensity * Si;
            }
        
            // Final color calculation
            finalColor = ambient + totalDiffuse + totalSpecular;
        }
        //clamp made the color to 0-1 value, but not seems neccessary some test cases with given color value
        //std::cout<<finalColor<<"final color before clamp is"<<std::endl;
        return finalColor.clamp();
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


  
    std::cout << "Rendering complete. Images saved as 'output_perspective.ppm' and 'output_parallel.ppm'." << std::endl;
    return 0;
}
