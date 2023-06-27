// Created by Tomoyuki Kawashima on 27/June/2023
// include guard
#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include <gmsh.h>
#include <vector3d.h>

#include <iostream>
#include <vector>

using namespace std;
namespace IM {
struct Node {
  size_t id;
  vec3d points[3];
};

struct Element {
  size_t id;
  size_t type;
  vector<int> nodeid;
  vec3d normal;
  double emissivity;
  double area;
};
}  // namespace IM
#endif  // _ELEMENT_H_
