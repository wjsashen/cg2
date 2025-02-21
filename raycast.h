#ifndef RAYCAST_H
#define RAYCAST_H
#include <fstream>
#include <vector>
#include <sstream>
#include <memory>
#include <cmath>
#include <limits>
#include <iostream>
#include <cmath>

#include <vector>
#include <string>
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
    float x, y, z;        // origin components
    float dx, dy, dz;     // direction components
    
    Ray(const Vec3d& o, const Vec3d& d)
        : x(o.x), y(o.y), z(o.z), 
          dx(d.x), dy(d.y), dz(d.z) {}
          
    Vec3d getOrigin() const { return Vec3d(x, y, z); }
    Vec3d getDirection() const { return Vec3d(dx, dy, dz); }
};


struct Light {
    Vec3d positionOrdir; //or is direction if it's not point light
    bool isPoint;
    float intensity;
    float c1, c2, c3;
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
    Vec3d center;
    double radius;
    MaterialColor material; 

    Sphere(Vec3d c, double r, MaterialColor mat)  // Accept MaterialColor
        : center(c), radius(r), material(mat) {}

    MaterialColor getColor() const override {
        return material;  // Return MaterialColor
    }

};

struct Cylinder : public Objects {
    Vec3d center;
    Vec3d dir;
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
//just attempts these 2
struct Intersection {
    bool hit;
    Vec3d point, normal;
    MaterialColor material;
};
class Plane : public Objects {
    public:
        Vec3d point;   
        Vec3d normal;   
        MaterialColor material;
    
        Plane(const Vec3d& p, const Vec3d& n, const MaterialColor& mat) 
            : point(p), normal(n.norm()), material(mat) {}
    
        MaterialColor getColor() const override {
            return material;
        }
    };
    #endif // RAYCAST_H