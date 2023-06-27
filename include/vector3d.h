// Created by Tomoyuki Kawashima on 27/June/2023
// include guard
#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include <cmath>
#include <vector>

namespace IM {
template <typename T>
class vec3d {
 public:
  T x;
  T y;
  T z;

  vec3d() : x(0), y(0), z(0) {}
  vec3d(const T _x, const T _y, const T _z) : x(_x), y(_y), z(_z) {}
  vec3d(const vec3d& obj) : x(obj.x), y(obj.y), z(obj.z) {}
  ~vec3d() {}
  const vec3d& operator=(const vec3d& obj) {
    x = obj.x;
    y = obj.y;
    z = obj.z;
    return *this;
  }
  const T& operator[](const int index) const {
    if (index == 0) {
      return x;
    } else if (index == 1) {
      return y;
    } else if (index == 2) {
      return z;
    } else {
      return x;
    }
  }
  T& operator[](const int index) {
    if (index == 0) {
      return x;
    } else if (index == 1) {
      return y;
    } else if (index == 2) {
      return z;
    } else {
      return x;
    }
  }
  vec3d operator+(const vec3d& obj) const {
    return vec3d(x + obj.x, y + obj.y, z + obj.z);
  }
  const vec3d& operator+=(const vec3d& obj) {
    x += obj.x;
    y += obj.y;
    z += obj.z;
    return *this;
  }
  vec3d operator+(const T& scal) { return vec3d(x + scal, y + scal, z + scal); }
  const vec3d& operator+=(const T& scal) {
    x += scal;
    y += scal;
    z += scal;
    return *this;
  }
  vec3d operator-(const vec3d& obj) {
    return vec3d(x - obj.x, y - obj.y, z - obj.z);
  }
  vec3d operator-() { return vec3d(-x, -y, -z); }
  vec3d operator-(const T& scal) { return vec3d(x - scal, y - scal, z - scal); }
  const vec3d& operator-=(const vec3d& obj) {
    x -= obj.x;
    y -= obj.y;
    z -= obj.z;
    return *this;
  }
  vec3d operator*(const vec3d& obj) {
    return vec3d(x * obj.x, y * obj.y, z * obj.z);
  }
  const vec3d& operator*=(const vec3d& obj) {
    x *= obj.x;
    y *= obj.y;
    z *= obj.z;
    return *this;
  }
  vec3d operator*(const T& scal) { return vec3d(x * scal, y * scal, z * scal); }
  const vec3d& operator*=(const T& scal) {
    x *= scal;
    y *= scal;
    z *= scal;
    return *this;
  }
  vec3d operator/(const vec3d& obj) {
    return vec3d(x / obj.x, y / obj.y, z / obj.z);
  }
  const vec3d& operator/=(const vec3d& obj) {
    x /= obj.x;
    y /= obj.y;
    z /= obj.z;
    return *this;
  }
  vec3d operator/(const T& scal) { return vec3d(x / scal, y / scal, z / scal); }
  const vec3d& operator/=(const T& scal) {
    x /= scal;
    y /= scal;
    z /= scal;
    return *this;
  }
  T dot(const vec3d& obj) { return x * obj.x + y * obj.y + z * obj.z; }
  vec3d cross(const vec3d& obj) {
    return vec3d(y * obj.z - z * obj.y, z * obj.x - x * obj.z,
                 x * obj.y - y * obj.x);
  }
};

// dot product
template <typename T>
T dot(const vec3d& obj1, const vec3d& obj2) {
  return obj1.x * obj2.x + obj1.y * obj2.y + obj1.z * obj2.z;
}
// cross product
vec3d cross(const vec3d& obj1, const vec3d& obj2) {
  return vec3d(obj1.y * obj2.z - obj1.z * obj2.y,
               obj1.z * obj2.x - obj1.x * obj2.z,
               obj1.x * obj2.y - obj1.y * obj2.x);
}
}  // namespace IM
#endif  // _VECTOR3D_H_
