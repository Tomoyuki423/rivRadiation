// Created by Tomoyuki Kawashima on 27/June/2023
// include guard
#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include <cmath>
#include <stdexcept>

namespace FemFy {

/**
 * @brief A template class for a 3D vector.
 *
 * @tparam T The type of the vector's coordinates.
 */
template <typename T>
class vec3d {
 public:
  T coord[3];

  /**
   * @brief Default constructor that initializes the vector to (0, 0, 0).
   */
  vec3d() : coord{0, 0, 0} {}

  /**
   * @brief Constructor that initializes the vector to (x, y, z).
   *
   * @param x The x-coordinate of the vector.
   * @param y The y-coordinate of the vector.
   * @param z The z-coordinate of the vector.
   */
  vec3d(T x, T y, T z) : coord{x, y, z} {}

  /**
   * @brief Copy constructor that initializes the vector to the same values as
   * another vector.
   *
   * @param obj The vector to copy.
   */
  vec3d(const vec3d& obj) : coord{obj.coord[0], obj.coord[1], obj.coord[2]} {}

  /**
   * @brief Assignment operator that sets the vector to the same values as
   * another vector.
   *
   * @param obj The vector to copy.
   * @return const vec3d& A reference to the modified vector.
   */
  const vec3d& operator=(const vec3d& obj) {
    coord[0] = obj.coord[0];
    coord[1] = obj.coord[1];
    coord[2] = obj.coord[2];
    return *this;
  }

  /**
   * @brief Subscript operator that returns the value of the coordinate at the
   * given index.
   *
   * @param index The index of the coordinate to retrieve.
   * @return const T& A reference to the coordinate value.
   * @throw std::out_of_range If the index is out of range.
   */
  const T& operator[](const int index) const {
    if (index < 0 || index >= 3) {
      throw std::out_of_range("Index out of range");
    }
    return coord[index];
  }

  /**
   * @brief Subscript operator that returns a reference to the coordinate at the
   * given index.
   *
   * @param index The index of the coordinate to retrieve.
   * @return T& A reference to the coordinate value.
   * @throw std::out_of_range If the index is out of range.
   */
  T& operator[](const int index) {
    if (index < 0 || index >= 3) {
      throw std::out_of_range("Index out of range");
    }
    return coord[index];
  }

  /**
   * @brief Addition operator that returns the sum of two vectors.
   *
   * @param obj The vector to add.
   * @return vec3d The sum of the two vectors.
   */
  vec3d operator+(const vec3d& obj) const {
    return vec3d(coord[0] + obj.coord[0], coord[1] + obj.coord[1],
                 coord[2] + obj.coord[2]);
  }

  /**
   * @brief Addition assignment operator that adds another vector to this
   * vector.
   *
   * @param obj The vector to add.
   * @return const vec3d& A reference to the modified vector.
   */
  const vec3d& operator+=(const vec3d& obj) {
    coord[0] += obj.coord[0];
    coord[1] += obj.coord[1];
    coord[2] += obj.coord[2];
    return *this;
  }

  /**
   * @brief Subtraction operator that returns the difference of two vectors.
   *
   * @param obj The vector to subtract.
   * @return vec3d The difference of the two vectors.
   */
  vec3d operator-(const vec3d& obj) const {
    return vec3d(coord[0] - obj.coord[0], coord[1] - obj.coord[1],
                 coord[2] - obj.coord[2]);
  }

  /**
   * @brief Subtraction assignment operator that subtracts another vector from
   * this vector.
   *
   * @param obj The vector to subtract.
   * @return const vec3d& A reference to the modified vector.
   */
  const vec3d& operator-=(const vec3d& obj) {
    coord[0] -= obj.coord[0];
    coord[1] -= obj.coord[1];
    coord[2] -= obj.coord[2];
    return *this;
  }

  /**
   * @brief Multiplication operator that returns the product of the vector and a
   * scalar.
   *
   * @param scal The scalar to multiply the vector by.
   * @return vec3d The product of the vector and the scalar.
   */
  vec3d operator*(const T& scal) const {
    return vec3d(coord[0] * scal, coord[1] * scal, coord[2] * scal);
  }

  /**
   * @brief Multiplication assignment operator that multiplies the vector by a
   * scalar.
   *
   * @param scal The scalar to multiply the vector by.
   * @return const vec3d& A reference to the modified vector.
   */
  const vec3d& operator*=(const T& scal) {
    coord[0] *= scal;
    coord[1] *= scal;
    coord[2] *= scal;
    return *this;
  }

  /**
   * @brief Division operator that returns the quotient of the vector and a
   * scalar.
   *
   * @param scal The scalar to divide the vector by.
   * @return vec3d The quotient of the vector and the scalar.
   */
  vec3d operator/(const T& scal) const {
    return vec3d(coord[0] / scal, coord[1] / scal, coord[2] / scal);
  }

  /**
   * @brief Division assignment operator that divides the vector by a scalar.
   *
   * @param scal The scalar to divide the vector by.
   * @return const vec3d& A reference to the modified vector.
   */
  const vec3d& operator/=(const T& scal) {
    coord[0] /= scal;
    coord[1] /= scal;
    coord[2] /= scal;
    return *this;
  }
};

/**
 * @brief Calculates the dot product of two vectors.
 *
 * @tparam T The type of the vector's coordinates.
 * @param obj1 The first vector.
 * @param obj2 The second vector.
 * @return T The dot product of the two vectors.
 */
template <typename T>
T dot(const vec3d<T>& obj1, const vec3d<T>& obj2) {
  return obj1.coord[0] * obj2.coord[0] + obj1.coord[1] * obj2.coord[1] +
         obj1.coord[2] * obj2.coord[2];
}

/**
 * @brief Calculates the cross product of two vectors.
 *
 * @tparam T The type of the vector's coordinates.
 * @param obj1 The first vector.
 * @param obj2 The second vector.
 * @return vec3d The cross product of the two vectors.
 */
template <typename T>
vec3d<T> cross(const vec3d<T>& obj1, const vec3d<T>& obj2) {
  return vec3d<T>(
      obj1.coord[1] * obj2.coord[2] - obj1.coord[2] * obj2.coord[1],
      obj1.coord[2] * obj2.coord[0] - obj1.coord[0] * obj2.coord[2],
      obj1.coord[0] * obj2.coord[1] - obj1.coord[1] * obj2.coord[0]);
}

}  // namespace FemFy

#endif  // _VECTOR3D_H_
