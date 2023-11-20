// Created by Tomoyuki Kawashima on 27/June/2023
// include guard
#ifndef _RADIATION_H_
#define _RADIATION_H_

#include <vector3d.h>

#include <iostream>
#include <vector>

using namespace std;
namespace FemFy {

// Define a function to interpolate the points on a surface
vec3d interpolate(const vec3d &points[3], double alpha, double beta,
                  double gamma) {
  vec3d interpolatedPoint;
  for (int i = 0; i < 3; ++i) {
    interpolatedPoint[i] = alpha * points[0].coord[i] +
                           beta * points[1].coord[i] +
                           gamma * points[2].coord[i];
  }
  return interpolatedPoint;
}

void generateInterpolatedPoints(const Point &A, const Point &B, const Point &C,
                                int n) {}
}

// / Define a function to check whether a ray intersects a surface
bool doesRayIntersectSurface(const vec3d &point1, const vec3d &point2,
                             const Surface &surface) {
  // Calculate the direction vector of the ray
  vec3d direction = point2 - point1;

  // Calculate the vectors of the triangle
  vec3d edge1 = surface.points[1] - surface.points[0]; // edge 1
  vec3d edge2 = surface.points[2] - surface.points[0]; // edge 2

  // Calculate the determinant
  vec3d h = cross(direction, edge2);
  double a = dot(edge1, h);

  // Check if the ray is parallel to the triangle
  if (a > -0.00001 && a < 0.00001) {
    return false;
  }

  // Calculate the inverse of a and check if the ray is pointing away from the
  // triangle (backface culling)
  double f = 1.0 / a;
  vec3d s = point1 - surface.points[0];

  // u parameter calculation and test if the intersection lies outside of the
  // edge 1 and edge 2 vectors of the triangle (the xy-plane) and if it does,
  // the ray misses the triangle and we can skip the rest of this test for
  // this triangle (continue) and move on to the next triangle.
  double u = f * dot(s, h);
  if (u < 0.0 || u > 1.0) {
    return false;
  }

  // v parameter calculation and test if the intersection lies outside of the
  // edge 1 and edge 2 vectors of the triangle (the xy-plane) and if it does,
  // the ray misses the triangle and we can skip the rest of this test for
  // this triangle (continue) and move on to the next triangle.
  vec3d q = cross(s, edge1);
  double v = f * dot(direction, q);
  if (v < 0.0 || u + v > 1.0) {
    return false;
  }

  // At this stage we can compute t to find out where the intersection point
  // is on the line
  double t = f * dot(edge2, q);
  if (t > 0.00001) { // ray intersection
    return true;
  }

  // This means that there is a line intersection but not a ray intersection
  return false;
}

// Define a function to check for obstructions between two points
bool isObstructed(const vec3d &point1, const vec3d &point2,
                  const std::vector<Surface> &allSurfaces) {
  // Loop over all the surfaces
  for (const Surface &surface : allSurfaces) {
    // Check whether the ray from point1 to point2 intersects this surface
    if (doesRayIntersectSurface(point1, point2, surface)) {
      // If it does, then there is an obstruction
      return true;
    }
  }

  // If we've checked all the surfaces and found no intersections, then there
  // is no obstruction
  return false;
}

// Define a function to calculate the view factor between two points
// F_12 = [(cosθ1 cosθ2) / (π r^2)] dA1 dA2
double
calculatePointToPointViewFactor(const vec3d &point1, const vec3d &normal1,
                                const vec3d &point2, const vec3d &normal2,
                                const std::vector<Surface> &allSurfaces) {
  if (isObstructed(point1, point2, allSurfaces)) {
    return 0.0;
  }
  // Calculate the distance between the points
  vec3d delta = point2 - point1;
  double r = dot(delta, delta);

  // Calculate the dot product of the normal at the first point and the
  // direction vector
  double dotProduct1 = dot(normal1, delta);

  // Calculate the dot product of the normal at the second point and the
  // direction vector
  double dotProduct2 = dot(normal2, delta);

  // Calculate the view factor
  double viewFactor = std::max(0.0, dotProduct1) * std::max(0.0, -dotProduct2) /
                      (M_PI * r * r * r * r);

  return viewFactor;
}

// Define a function to calculate the view factor between two surfaces
double calculateViewFactor(const Surface &surface1, const Surface &surface2,
                           const std::vector<Surface> &allSurfaces) {
  // Number of points to use for the numerical integration
  int numPoints = 5;

  // Initialize the view factor to zero
  double viewFactor = 0.0;

  // Loop over each point on surface 1
  for (int i = 0; i <= numPoints; ++i) {
    for (int j = 0; j <= numPoints - i; ++j) {
      int k = numPoints - i - j;
      double alpha = (double)i / numPoints;
      double beta = (double)j / numPoints;
      double gamma = (double)k / numPoints;

      // Calculate the position and normal of this point
      vec3d point1 = interpolate(surface1.points, alpha, beta, gamma);
      vec3d normal1 = surface1.normal;

      // Loop over each point on surface 2
      for (int l = 0; l <= numPoints; ++l) {
        for (int m = 0; m <= numPoints - l; ++m) {
          int n = numPoints - l - m;
          alpha = (double)l / numPoints;
          beta = (double)m / numPoints;
          gamma = (double)n / numPoints;
          // Calculate the position and normal of this point
          vec3d point1 = interpolate(surface2.points, alpha, beta, gamma);
          vec3d normal2 = surface2.normal;

          // Add the contribution of this pair of points to the view factor
          viewFactor += calculatePointToPointViewFactor(point1, normal1, point2,
                                                        normal2, allSurfaces) /
                        (numPoints * numPoints);
        }
      }
    }
  }

  return viewFactor;
}

// Define a function to calculate the radiative heat flux between two surfaces
double calculateRadiativeHeatFlux(const Surface &surface1,
                                  const Surface &surface2) {
  double viewFactor = calculateViewFactor(surface1, surface2);
  double radiativeHeatFlux =
      surface1.emissivity * surface1.area * viewFactor *
      (std::pow(surface1.temperature, 4) - std::pow(surface2.temperature, 4));
  return radiativeHeatFlux;
}

// Define a function to calculate the total radiative heat transfer from one
// surface to all others
double
calculateTotalRadiativeHeatTransfer(const Surface &surface,
                                    const std::vector<Surface> &allSurfaces) {
  double totalRadiativeHeatTransfer = 0.0;
  for (const Surface &otherSurface : allSurfaces) {
    if (&surface != &otherSurface) {
      totalRadiativeHeatTransfer +=
          calculateRadiativeHeatFlux(surface, otherSurface);
    }
  }
  return totalRadiativeHeatTransfer;
}
} // namespace FemFy
#endif // _RADIATION_H_
