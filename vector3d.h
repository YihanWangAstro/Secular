#ifndef GENVECT_H
#define GENVECT_H

#include <math.h>
#include <iostream>

    struct Vec3d {
    public:
        /* Typedef */
        using value_type = double;
        /* Typedef */

        value_type x{0};
        value_type y{0};
        value_type z{0};

        Vec3d() = default;

        Vec3d(value_type vx, value_type vy, value_type vz) : x(vx), y(vy), z(vz) {}

        Vec3d(const Vec3d &v) = default;

        Vec3d(Vec3d&&) = default;

        Vec3d operator=(Vec3d const& v){
            x = v.x, y = v.y, z = v.z; 
            return *this;
        }

        Vec3d operator=(Vec3d && v){
            x = v.x, y = v.y, z = v.z; 
            return *this;
        }

        /** @brief Addition by wise */
        inline Vec3d operator+(const Vec3d &v) const {
            return Vec3d(x + v.x, y + v.y, z + v.z);
        }

        /** @brief Subtraction by wise */
        inline Vec3d operator-(const Vec3d &v) const {
            return Vec3d(x - v.x, y - v.y, z - v.z);
        }

        /** @brief Product by wise */
        inline Vec3d operator*(const Vec3d &v) const {
            return Vec3d(x * v.x, y * v.y, z * v.z);
        }

        /** @brief Divition by wise */
        inline Vec3d operator/(const Vec3d &v) const {
            return Vec3d(x / v.x, y / v.y, z / v.z);
        }

        /** @brief Add scalar by wise */
        inline Vec3d operator+(const value_type c) const {
            return Vec3d(x + c, y + c, z + c);
        }

        /** @brief Subtract scalar by wise */
        inline Vec3d operator-(const value_type c) const {
            return Vec3d(x - c, y - c, z - c);
        }

        /** @brief Multiply scalar by wise */
        inline Vec3d operator*(const value_type c) const {
            return Vec3d(x * c, y * c, z * c);
        }

        /** @brief Divide scalar by wise */
        inline Vec3d operator/(const value_type c) const {
            return Vec3d(x / c, y / c, z / c);
        }

        /** @brief Opposite vector */
        inline Vec3d operator-() const {
            return Vec3d(-x, -y, -z);
        }

        inline const Vec3d &operator+=(const Vec3d &v) {
            x += v.x, y += v.y, z += v.z;
            return *this;
        }

        inline const Vec3d &operator-=(const Vec3d &v) {
            x -= v.x, y -= v.y, z -= v.z;
            return *this;
        }

        inline const Vec3d &operator*=(const Vec3d &v) {
            x *= v.x, y *= v.y, z *= v.z;
            return *this;
        }

        inline const Vec3d &operator/=(const Vec3d &v) {
            x /= v.x, y /= v.y, z /= v.z;
            return *this;
        }

        inline const Vec3d &operator+=(const value_type c) {
            x += c, y += c, z += c;
            return *this;
        }

        inline const Vec3d &operator-=(const value_type c) {
            x -= c, y -= c, z -= c;
            return *this;
        }

        inline const Vec3d &operator*=(const value_type c) {
            x *= c, y *= c, z *= c;
            return *this;
        }

        inline const Vec3d &operator/=(const value_type c) {
            x /= c, y /= c, z /= c;
            return *this;
        }

        friend Vec3d operator+(const value_type c, const Vec3d &v) {
            return Vec3d(v.x + c, v.y + c, v.z + c);
        }

        friend Vec3d operator-(const value_type c, const Vec3d &v) {
            return Vec3d(c - v.x, c - v.y, c - v.z);
        }

        friend Vec3d operator*(const value_type c, const Vec3d &v) {
            return Vec3d(v.x * c, v.y * c, v.z * c);
        }

        friend Vec3d operator/(const value_type c, const Vec3d &v) {
            return Vec3d(c / v.x, c / v.y, c / v.z);
        }

        /** @brief Output to ostream */
        friend std::ostream &operator<<(std::ostream &output, const Vec3d &v) {
            output << v.x << ' ' << v.y << ' ' << v.z;
            return output;
        }

        /** @brief Input from istream */
        friend std::istream &operator>>(std::istream &input, Vec3d &v) {
            input >> v.x >> v.y >> v.z;
            return input;
        }
    };


/** @brief Calculate the inner product of two vectors */

    inline double dot(const Vec3d &v1, const Vec3d &v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    inline Vec3d unit(const Vec3d& v) {
        double s = 1.0/sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        return Vec3d(v.x * s, v.y * s, v.z*s);
    }
   
    inline double norm(const Vec3d &v) {
        return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    
    inline double re_norm(const Vec3d &v) {
        return 1.0/sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    
    inline double norm2(const Vec3d &v) {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }


/** @brief Calculate the cross product of two vectors */
    
    inline Vec3d cross(const Vec3d &v1, const Vec3d &v2) {
        return Vec3d(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }


#endif