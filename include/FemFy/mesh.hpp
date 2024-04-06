// Created by Tomoyuki Kawashima on 27/June/2023
// include guard
#ifndef _MESH_H_
#define _MESH_H_

// #include <vcruntime.h>

#include <FemFy/FemFy.hpp>
#include <FemFy/utils.hpp>

namespace FemFy {

enum class GmshElementType {
  POINT_1NODE = 15,
  // 1D Elements
  LINE_2NODE = 1,
  LINE_3NODE = 8,

  // 2D Elements
  TRIANGLE_3NODE = 2,
  TRIANGLE_6NODE = 9,
  QUAD_4NODE = 3,
  QUAD_8NODE = 16,
  QUAD_9NODE = 10,

  // 3D Elements
  TETRAHEDRON_4NODE = 4,
  TETRAHEDRON_10NODE = 11,
  HEXAHEDRON_8NODE = 5,
  HEXAHEDRON_20NODE = 12,
  HEXAHEDRON_27NODE = 29,
  PRISM_6NODE = 6,
  PRISM_15NODE = 13,
  PRISM_18NODE = 21,
  PYRAMID_5NODE = 7,
  PYRAMID_13NODE = 14,
  PYRAMID_14NODE = 32,

  // Other types can be added as needed
};

template <typename T>
class Node {
 protected:
  size_t id;
  vec3d<T> pos;

 public:
  Node() : id(0), pos(0.0, 0.0, 0.0){};  // Initialize id in default constructor
  explicit Node(size_t _id, T x, T y, T z) : id(_id) {
    pos.coord[0] = x;
    pos.coord[1] = y;
    pos.coord[2] = z;
  }  // Added constructor with parameters
  explicit Node(size_t _id, const vec3d<T>& p)
      : id(_id), pos(p) {}  // Added constructor with parameters
  ~Node() = default;
  vec3d<T> getPos() const { return pos; }
  T getPos(int i) const {
    if (i > 3) {
      std::runtime_error("Here is getPos function. i is larger than 3");
    }
    return pos.coord[i];
  }
  void setPos(const vec3d<T>& p) { pos = p; }
  size_t getId() const { return id; }  // Getter for id
  //   std::vector<T> phi;
  //   std::map<std::string, int> variable;
  double temperature;
};

class Element {
 protected:
  size_t id;
  GmshElementType type;  // Changed type to GmshElementType enum class
  std::vector<size_t> nodeIds;
  const int expectedNodeCount;
  double shapeFuctions;
  double shapeFuctionsDerivatives;

 public:
  GmshElementType getType() const { return type; }
  Element(size_t id, GmshElementType type, const std::vector<size_t>& nodeIds,
          int expectedNodeCount)
      : id(id),
        type(type),
        nodeIds(nodeIds),
        expectedNodeCount(expectedNodeCount) {}
  virtual ~Element() = default;

  virtual size_t getNodeCount() const = 0;
  // virtual size_t getNodeIDs(const int i) const = 0;
  int getNodeIds(const int i) const {
    if (i > expectedNodeCount) {
      std::runtime_error("i is larger than expectedNodeCount");
    }
    return nodeIds[i];
  }
  // 形状関数の計算
  virtual std::vector<double> computeShapeFunctions(
      const double& xi, const double& eta = 0,
      const double& zeta = 0) const = 0;

  // 形状関数の導関数の計算
  virtual std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& xi, const double& eta = 0,
      const double& zeta = 0) const = 0;

  virtual void computeDimensionalProperties(
      std::map<int, Node<double>> nodes) = 0;

  size_t getId() const { return id; }
  virtual double getDimensionalProperties() const = 0;

  virtual void printType() const = 0;

  virtual void printInfo() const {
    std::cout << "Element ID: " << id << ", Type: " << static_cast<int>(type)
              << ", ";
    printType();
    std::cout << ", Node IDs: ";
    for (size_t nodeId : nodeIds) {
      std::cout << nodeId << " ";
    }
    std::cout << std::endl;
  }
};
class Point : public Element {
 public:
  Point(size_t id, const std::vector<size_t>& nodeIds)
      : Element(id, GmshElementType::POINT_1NODE, nodeIds, 1) {}

  void printType() const override { std::cout << "Point"; }
  size_t getNodeCount() const override { return expectedNodeCount; }
  double getDimensionalProperties() const {
    std::cout << "sorry not implemented yet about dimension properties in Point"
              << std::endl;
  };
  // 形状関数の計算
  void computeDimensionalProperties(std::map<int, Node<double>> nodes) {
    std::cout << "sorry not implemented yet about dimension properties in Point"
              << std::endl;
  }
  std::vector<double> computeShapeFunctions(
      const double& xi, const double& eta = 0,
      const double& zeta = 0) const override {
    return {1.0};
  }

  // 形状関数の導関数の計算
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& xi, const double& eta = 0,
      const double& zeta = 0) const override {
    return {{0}, {0}};
  }
};

class Line : public Element {
 public:
  double length;
  Line(size_t id, GmshElementType type, const std::vector<size_t>& nodeIds,
       int expectedNodeCount)
      : Element(id, type, nodeIds, expectedNodeCount) {}
  virtual void printType() const override { std::cout << "Line"; }
  double getDimensionalProperties() const { return getLength(); }
  double getLength() const { return length; }
  // 1次元要素特有のプロパティやメソッド
  void computeDimensionalProperties(std::map<int, Node<double>> nodes) {
    std::cout << "sorry not implemented yet about dimension properties in Line"
              << std::endl;
  }
};

class Surface : public Element {
 public:
  double area;
  Surface(size_t id, GmshElementType type, const std::vector<size_t>& nodeIds,
          int expectedNodeCount)
      : Element(id, type, nodeIds, expectedNodeCount) {}
  virtual void printType() const override { std::cout << "Surface"; }
  double getDimensionalProperties() const { return getArea(); }
  double getArea() const { return area; }

  // 2次元要素特有のプロパティやメソッド
  void computeDimensionalProperties(std::map<int, Node<double>> nodes) {
    computeArea(nodes);
  };
  virtual void computeArea(std::map<int, Node<double>> nodes) = 0;
};

class Volume : public Element {
 public:
  double volume;
  Volume(size_t id, GmshElementType type, const std::vector<size_t>& nodeIds,
         int expectedNodeCount)
      : Element(id, type, nodeIds, expectedNodeCount) {}
  virtual void printType() const override { std::cout << "Volume"; }
  double getDimensionalProperties() const { return getVolume(); }
  double getVolume() const { return volume; }
  // 3次元要素特有のプロパティやメソッド
  void computeDimensionalProperties(std::map<int, Node<double>> nodes) {
    computeVolume(nodes);
  }
  virtual void computeVolume(std::map<int, Node<double>> nodes) = 0;
};

// 線要素
class LineElement : public Line {
 public:
  // static constexpr size_t expectedNodeCount = 2;

  LineElement(size_t id, const std::vector<size_t>& nodeIds)
      : Line(id, GmshElementType::LINE_2NODE, nodeIds, 2) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument("LineElement requires exactly 2 nodes");
    }
  }

  size_t getNodeCount() const override { return expectedNodeCount; }

  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(1);
    N = {0.5 * (1 - r), 0.5 * (1 + r)};
    return N;
  }
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(1, std::vector<double>(2));
    dNdrst[0] = {-0.5, 0.5};
    return dNdrst;
  }
};

class QuadraticLineElement : public Line {
 public:
  // static constexpr size_t expectedNodeCount = 3;
  QuadraticLineElement(size_t id, const std::vector<size_t>& nodeIds)
      : Line(id, GmshElementType::LINE_3NODE, nodeIds, 3) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument(
          "QuadraticLineElement requires exactly 3 nodes");
    }
  }

  size_t getNodeCount() const override { return expectedNodeCount; }
  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(3);
    N[0] = 0.5 * r * (r - 1);
    N[1] = 1 - r * r;
    N[2] = 0.5 * r * (r + 1);
    return N;
  }
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdr(1, std::vector<double>(3));
    dNdr[0] = {r - 0.5, -2 * r, r + 0.5};
    return dNdr;
  }
};

// 三角形要素
class TriangleElement : public Surface {
 public:
  // static constexpr size_t expectedNodeCount = 3;

  TriangleElement(size_t id, const std::vector<size_t>& nodeIds)
      : Surface(id, GmshElementType::TRIANGLE_3NODE, nodeIds, 3) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument("TriangleElement requires exactly 3 nodes");
    }
  }

  size_t getNodeCount() const override { return expectedNodeCount; }

  void computeArea(std::map<int, Node<double>> nodes) override {
    vec3d<double> v1 = nodes[nodeIds[1]].getPos() - nodes[nodeIds[0]].getPos();
    vec3d<double> v2 = nodes[nodeIds[2]].getPos() - nodes[nodeIds[0]].getPos();
    area = 0.5 * v1.cross(v2).norm();
    if (area < 0) {
      std::cout << "Area is negative" << std::endl;
    }
  }

  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(3);
    N[0] = 1 - r - s;
    N[1] = r;
    N[2] = s;
    return N;
  }  // 形状関数の計算 (r, s, t) = (ξ, η, ζ) に対して  N1 = 1 - ξ - η, N2 = ξ,
     // N3 = η を返す  (三角形要素)
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(3, std::vector<double>(2));
    dNdrst[0] = {-1, -1};
    dNdrst[1] = {1, 0};
    dNdrst[2] = {0, 1};
    return dNdrst;
  }  // 形状関数の導関数の計算 (r, s, t) = (ξ, η, ζ) に対して  dN1/dξ = -1,
     // dN1/dη = -1, dN2/dξ = 1, dN2/dη = 0, dN3/dξ = 0, dN3/dη = 1 を返す
     // (三角形要素)
};

// 四角形要素
class QuadrangleElement : public Surface {
 public:
  // static constexpr size_t expectedNodeCount = 4;

  QuadrangleElement(size_t id, const std::vector<size_t>& nodeIds)
      : Surface(id, GmshElementType::QUAD_4NODE, nodeIds, 4) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument("QuadrangleElement requires exactly 4 nodes");
    }
  }
  size_t getNodeCount() const override { return expectedNodeCount; }

  void computeArea(std::map<int, Node<double>> nodes) override {
    // vec3d<double> v1 = nodes[nodeIds[1]].getPos() -
    // nodes[nodeIds[0]].getPos(); vec3d<double> v2 =
    // nodes[nodeIds[2]].getPos()
    // - nodes[nodeIds[0]].getPos(); area = 0.5 * v1.cross(v2).norm();
    throw std::runtime_error(
        "Not implemented yet about compute area in QuadrangleElement");
  }

  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(4);
    N[0] = 0.25 * (1 - r) * (1 - s);
    N[1] = 0.25 * (1 + r) * (1 - s);
    N[2] = 0.25 * (1 + r) * (1 + s);
    N[3] = 0.25 * (1 - r) * (1 + s);
    return N;
  }
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(4, std::vector<double>(2));
    dNdrst[0] = {-0.25 * (1 - s), -0.25 * (1 - r)};
    dNdrst[1] = {0.25 * (1 - s), -0.25 * (1 + r)};
    dNdrst[2] = {0.25 * (1 + s), 0.25 * (1 + r)};
    dNdrst[3] = {-0.25 * (1 + s), 0.25 * (1 - r)};
    return dNdrst;
  }
};

// 二次三角形要素
class QuadraticTriangleElement : public Surface {
 public:
  // static constexpr size_t expectedNodeCount = 6;

  QuadraticTriangleElement(size_t id, const std::vector<size_t>& nodeIds)
      : Surface(id, GmshElementType::TRIANGLE_6NODE, nodeIds, 6) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument(
          "QuadraticTriangleElement requires exactly 6 nodes");
    }
  }

  size_t getNodeCount() const override { return expectedNodeCount; }

  void computeArea(std::map<int, Node<double>> nodes) override {
    vec3d<double> v1 = nodes[nodeIds[1]].getPos() - nodes[nodeIds[0]].getPos();
    vec3d<double> v2 = nodes[nodeIds[2]].getPos() - nodes[nodeIds[0]].getPos();
    area = 0.5 * v1.cross(v2).norm();
    // throw std::runtime_error(
    //     "Not implemented yet about compute area in "
    //     "QuandraticTriangleElement");
  }
  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(6);
    N[0] = (1 - r - s) * (1 - 2 * r - 2 * s);
    N[1] = r * (2 * r - 1);
    N[2] = s * (2 * s - 1);
    N[3] = 4 * r * (1 - r - s);
    N[4] = 4 * r * s;
    N[5] = 4 * s * (1 - r - s);
    return N;
  }  // 形状関数の計算 (r, s, t) = (ξ, η, ζ) に対して  N1 = (1 - ξ - η)(1 - 2ξ
     // - 2η),  N2 = ξ(2ξ - 1),  N3 = η(2η - 1),  N4 = 4ξ(1 - ξ - η),  N5 =
     // 4ξη, N6 = 4η(1 - ξ - η) を返す  (二次三角形要素)
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(6, std::vector<double>(2));
    dNdrst[0] = {-3 + 4 * r + 4 * s, -3 + 4 * r + 4 * s};
    dNdrst[1] = {4 * r - 1, 0};
    dNdrst[2] = {0, 4 * s - 1};
    dNdrst[3] = {4 - 8 * r - 4 * s, -4 * r};
    dNdrst[4] = {4 * s, 4 * r};
    dNdrst[5] = {-4 * s, 4 - 4 * r - 8 * s};
    return dNdrst;
  }  // 形状関数の導関数の計算 (r, s, t) = (ξ, η, ζ) に対して  dN1/dξ = -3 +
     // 4ξ
     // + 4η,  dN1/dη = -3 + 4ξ + 4η,  dN2/dξ = 4ξ - 1,  dN2/dη = 0,  dN3/dξ =
     // 0,  dN3/dη = 4η - 1,  dN4/dξ = 4 - 8ξ - 4η,  dN4/dη = -4ξ,  dN5/dξ =
     // 4η, dN5/dη = 4ξ,  dN6/dξ = -4η,  dN6/dη = 4 - 4ξ - 8η を返す
     // (二次三角形要素)
};

// 二次四角形要素
class QuadraticQuadrangleElement : public Surface {
 public:
  // static constexpr size_t expectedNodeCount = 8;

  QuadraticQuadrangleElement(size_t id, const std::vector<size_t>& nodeIds)
      : Surface(id, GmshElementType::QUAD_8NODE, nodeIds, 8) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument(
          "QuadraticQuadrangleElement requires exactly 8 nodes");
    }
  }

  size_t getNodeCount() const override { return expectedNodeCount; }

  void computeArea(std::map<int, Node<double>> nodes) override {
    // vec3d<double> v1 = nodes[nodeIds[1]].getPos() -
    // nodes[nodeIds[0]].getPos(); vec3d<double> v2 =
    // nodes[nodeIds[2]].getPos()
    // - nodes[nodeIds[0]].getPos(); area = 0.5 * v1.cross(v2).norm();
    throw std::runtime_error(
        "Not implemented yet about compute area in "
        "QuandraticQuadraticQuadrangleElement");
  }

  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(8);
    N[0] = 0.25 * (1 - r) * (1 - s) * (-r - s - 1);
    N[1] = 0.25 * (1 + r) * (1 - s) * (r - s - 1);
    N[2] = 0.25 * (1 + r) * (1 + s) * (r + s - 1);
    N[3] = 0.25 * (1 - r) * (1 + s) * (-r + s - 1);
    N[4] = 0.5 * (1 - r * r) * (1 - s);
    N[5] = 0.5 * (1 + r) * (1 - s * s);
    N[6] = 0.5 * (1 - r * r) * (1 + s);
    N[7] = 0.5 * (1 - r) * (1 - s * s);
    return N;
  }  // 形状関数の計算 (r, s, t) = (ξ, η, ζ) に対して  N1 = 0.25(1 - ξ)(1 -
     // η)(-ξ
     // - η - 1),  N2 = 0.25(1 + ξ)(1 - η)(ξ - η - 1),  N3 = 0.25(1 + ξ)(1 +
     // η)(ξ + η - 1),  N4 = 0.25(1 - ξ)(1 + η)(-ξ + η - 1),  N5 = 0.5(1 -
     // ξ^2)(1
     // - η),  N6 = 0.5(1 + ξ)(1 - η^2),  N7 = 0.5(1 - ξ^2)(1 + η),  N8 =
     // 0.5(1
     // - ξ)(1 - η^2) を返す  (二次四角形要素)
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(8, std::vector<double>(2));
    dNdrst[0] = {0.25 * (s - 1) * (2 * r + s), 0.25 * (r - 1) * (r + 2 * s)};
    dNdrst[1] = {0.25 * (s - 1) * (2 * r - s), 0.25 * (r + 1) * (2 * r - s)};
    dNdrst[2] = {0.25 * (s + 1) * (2 * r + s), 0.25 * (r + 1) * (2 * r + s)};
    dNdrst[3] = {0.25 * (s + 1) * (2 * r - s), 0.25 * (r - 1) * (2 * r - s)};
    dNdrst[4] = {-r * (1 - s), 0.5 * (1 - r * r)};
    dNdrst[5] = {0.5 * (1 - s * s), -s * (1 + r)};
    dNdrst[6] = {-r * (1 + s), 0.5 * (1 - r * r)};
    dNdrst[7] = {0.5 * (1 - s * s), -s * (1 - r)};
    return dNdrst;
  }  // 形状関数の導関数の計算 (r, s, t) = (ξ, η, ζ) に対して  dN1/dξ =
     // 0.25(η - 1)(2ξ + η),  dN1/dη = 0.25(ξ - 1)(ξ + 2η),  dN2/dξ =
     // 0.25(η - 1)(2ξ - η),  dN2/dη = 0.25(ξ + 1)(2ξ - η),  dN3/dξ =
     // 0.25(η + 1)(2ξ + η),  dN3/dη = 0.25(ξ + 1)(2ξ + η),  dN4/dξ =
     // 0.25(η + 1)(2ξ - η),  dN4/dη = 0.25(ξ - 1)(2ξ - η),  dN5/dξ = -r(1 -
     // η), dN5/dη = 0.5(1 - r^2),  dN6/dξ = 0.5(1 - η^2),  dN6/dη = -s(1 +
     // r), dN7/dξ = -r(1 + η),  dN7/dη = 0.5(1 - r^2),  dN8/dξ = 0.5(1 -
     // η^2), dN8/dη = -s(1 - r) を返す  (二次四角形要素)
};

// 四面体要素
class TetrahedronElement : public Volume {
 public:
  // static constexpr size_t expectedNodeCount = 4;

  TetrahedronElement(size_t id, const std::vector<size_t>& nodeIds)
      : Volume(id, GmshElementType::TETRAHEDRON_4NODE, nodeIds, 4) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument(
          "TetrahedronElement requires exactly 4 nodes");
    }
  }

  size_t getNodeCount() const override { return expectedNodeCount; }
  void computeVolume(std::map<int, Node<double>> nodes) override {
    std::vector<std::vector<double>> matrix(3, std::vector<double>(3));
    vec3d<double> v1 = nodes[nodeIds[1]].getPos() - nodes[nodeIds[0]].getPos();
    vec3d<double> v2 = nodes[nodeIds[2]].getPos() - nodes[nodeIds[0]].getPos();
    vec3d<double> v3 = nodes[nodeIds[3]].getPos() - nodes[nodeIds[0]].getPos();
    for (int i = 0; i < 3; i++) {
      matrix[0][i] = v1.coord[i];
      matrix[1][i] = v2.coord[i];
      matrix[2][i] = v3.coord[i];
    }
    volume = 1.0 / 6.0 * calculateDeterminant3x3(matrix);
    if (volume < 0) {
      std::cout << "Area is negative" << std::endl;
    }
  }
  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(4);
    N[0] = 1 - r - s - t;
    N[1] = r;
    N[2] = s;
    N[3] = t;
    return N;
  }  // 形状関数の計算 (r, s, t) = (ξ, η, ζ) に対して  N1 = 1 - ξ - η - ζ, N2
     // = ξ, N3 = η, N4 = ζ を返す  (四面体要素)
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(4, std::vector<double>(3));
    dNdrst[0] = {-1, -1, -1};
    dNdrst[1] = {1, 0, 0};
    dNdrst[2] = {0, 1, 0};
    dNdrst[3] = {0, 0, 1};
    return dNdrst;
  }  // 形状関数の導関数の計算 (r, s, t) = (ξ, η, ζ) に対して  dN1/dξ = -1,
     // dN1/dη = -1, dN1/dζ = -1, dN2/dξ = 1, dN2/dη = 0, dN2/dζ = 0, dN3/dξ =
     // 0, dN3/dη = 1, dN3/dζ = 0, dN4/dξ = 0, dN4/dη = 0, dN4/dζ = 1 を返す
     // (四面体要素)
};

// 二次テトラヘドロン要素
class QuadraticTetrahedronElement : public Volume {
 public:
  QuadraticTetrahedronElement(size_t id, const std::vector<size_t>& nodeIds)
      : Volume(id, GmshElementType::TETRAHEDRON_10NODE, nodeIds, 10) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument(
          "QuadraticTetrahedronElement requires exactly 10 nodes");
    }
  }
  size_t getNodeCount() const override { return expectedNodeCount; }
  void computeVolume(std::map<int, Node<double>> nodes) override {
    std::vector<std::vector<double>> matrix(3, std::vector<double>(3));
    vec3d<double> v1 = nodes[nodeIds[1]].getPos() - nodes[nodeIds[0]].getPos();
    vec3d<double> v2 = nodes[nodeIds[2]].getPos() - nodes[nodeIds[0]].getPos();
    vec3d<double> v3 = nodes[nodeIds[3]].getPos() - nodes[nodeIds[0]].getPos();
    for (int i = 0; i < 3; i++) {
      matrix[0][i] = v1.coord[i];
      matrix[1][i] = v2.coord[i];
      matrix[2][i] = v3.coord[i];
    }
    volume = 1.0 / 6.0 * calculateDeterminant3x3(matrix);
    if (volume < 0) {
      std::cout << "Area is negative" << std::endl;
    }
    // throw std::runtime_error(
    //     "Not implemented yet about compute Volume in "
    //     "QuandraticTeraHedronElement");
  }
  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(10);
    N[0] = (1 - r - s - t) * (2 * (1 - r - s - t) - 1);
    N[1] = r * (2 * r - 1);
    N[2] = s * (2 * s - 1);
    N[3] = t * (2 * t - 1);
    N[4] = 4 * r * (1 - r - s - t);
    N[5] = 4 * r * s;
    N[6] = 4 * s * (1 - r - s - t);
    N[7] = 4 * t * (1 - r - s - t);
    N[8] = 4 * s * t;
    N[9] = 4 * r * t;
    return N;
  }

  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(10, std::vector<double>(3));
    // dN/dx derivatives
    dNdrst[0] = {4 * r + 4 * s + 4 * t - 3, 4 * r + 4 * s + 4 * t - 3,
                 4 * r + 4 * s + 4 * t - 3};
    dNdrst[1] = {4 * r - 1, 0, 0};
    dNdrst[2] = {0, 4 * s - 1, 0};
    dNdrst[3] = {0, 0, 4 * t - 1};
    dNdrst[4] = {4 - 8 * r - 4 * s - 4 * t, -4 * r, -4 * r};
    dNdrst[5] = {4 * s, 4 * r, 0};
    dNdrst[6] = {-4 * s, 4 - 4 * r - 8 * s - 4 * t, -4 * s};
    dNdrst[7] = {-4 * t, -4 * t, 4 - 4 * r - 4 * s - 8 * t};
    dNdrst[8] = {0, 4 * t, 4 * s};
    dNdrst[9] = {4 * t, 0, 4 * r};
    return dNdrst;
  }
};
// プリズム要素
class PrismElement : public Volume {
 public:
  // static constexpr size_t expectedNodeCount = 6;

  PrismElement(size_t id, const std::vector<size_t>& nodeIds)
      : Volume(id, GmshElementType::PRISM_6NODE, nodeIds, 6) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument("PrismElement requires exactly 6 nodes");
    }
  }

  size_t getNodeCount() const override { return expectedNodeCount; }
  void computeVolume(std::map<int, Node<double>> nodes) override {
    //   vec3d<double> v1 = nodes[nodeIds[1]].getPos() -
    //   nodes[nodeIds[0]].getPos(); vec3d<double> v2 =
    //   nodes[nodeIds[2]].getPos() - nodes[nodeIds[0]].getPos(); area = 0.5 *
    //   v1.cross(v2).norm();
    throw std::runtime_error(
        "Not implemented yet about compute Volume in "
        "PrismElement");
  }
  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(6);
    // 形状関数の定義
    N[0] = (1 - r - s) * (1 - t);
    N[1] = r * (1 - t);
    N[2] = s * (1 - t);
    N[3] = (1 - r - s) * (1 + t);
    N[4] = r * (1 + t);
    N[5] = s * (1 + t);

    return N;
  }

  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(6, std::vector<double>(3));
    // 形状関数の導関数の定義
    dNdrst[0] = {-(1 - t), -(1 - t), -(1 - r - s)};
    dNdrst[1] = {(1 - t), 0, -r};
    dNdrst[2] = {0, (1 - t), -s};
    dNdrst[3] = {-(1 + t), -(1 + t), 1 - r - s};
    dNdrst[4] = {(1 + t), 0, r};
    dNdrst[5] = {0, (1 + t), s};

    return dNdrst;
  }
};

// // 二次プリズム要素
// class QuadraticPrismElement : public Volume {
//  public:
//   static constexpr size_t expectedNodeCount = 15;

//   QuadraticPrismElement(size_t id, const std::vector<size_t>& nodeIds)
//       : Volume(id, GmshElementType::PRISM_15NODE, nodeIds) {
//     if (nodeIds.size() != expectedNodeCount) {
//       throw std::invalid_argument(
//           "QuadraticPrismElement requires exactly 15 nodes");
//     }
//   }

//   size_t getNodeCount() const override { return expectedNodeCount; }
// };

// ヘキサヘドロン要素
class HexahedronElement : public Volume {
 public:
  // static constexpr size_t expectedNodeCount = 8;

  HexahedronElement(size_t id, const std::vector<size_t>& nodeIds)
      : Volume(id, GmshElementType::HEXAHEDRON_8NODE, nodeIds, 8) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument("HexahedronElement requires exactly 8 nodes");
    }
  }
  void computeVolume(std::map<int, Node<double>> nodes) override {
    //   vec3d<double> v1 = nodes[nodeIds[1]].getPos() -
    //   nodes[nodeIds[0]].getPos(); vec3d<double> v2 =
    //   nodes[nodeIds[2]].getPos() - nodes[nodeIds[0]].getPos(); area = 0.5 *
    //   v1.cross(v2).norm();
    throw std::runtime_error(
        "Not implemented yet about compute Volume in "
        "HexahedronElement");
  }
  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(8);
    N[0] = 0.125 * (1 - r) * (1 - s) * (1 - t);
    N[1] = 0.125 * (1 + r) * (1 - s) * (1 - t);
    N[2] = 0.125 * (1 + r) * (1 + s) * (1 - t);
    N[3] = 0.125 * (1 - r) * (1 + s) * (1 - t);
    N[4] = 0.125 * (1 - r) * (1 - s) * (1 + t);
    N[5] = 0.125 * (1 + r) * (1 - s) * (1 + t);
    N[6] = 0.125 * (1 + r) * (1 + s) * (1 + t);
    N[7] = 0.125 * (1 - r) * (1 + s) * (1 + t);
    return N;
  }

  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(8, std::vector<double>(3));
    dNdrst[0] = {-0.125 * (1 - s) * (1 - t), -0.125 * (1 - r) * (1 - t),
                 -0.125 * (1 - r) * (1 - s)};
    dNdrst[1] = {0.125 * (1 - s) * (1 - t), -0.125 * (1 + r) * (1 - t),
                 -0.125 * (1 + r) * (1 - s)};
    dNdrst[2] = {0.125 * (1 + s) * (1 - t), 0.125 * (1 + r) * (1 - t),
                 -0.125 * (1 + r) * (1 + s)};
    dNdrst[3] = {-0.125 * (1 + s) * (1 - t), 0.125 * (1 - r) * (1 - t),
                 -0.125 * (1 - r) * (1 + s)};
    dNdrst[4] = {-0.125 * (1 - s) * (1 + t), -0.125 * (1 - r) * (1 + t),
                 0.125 * (1 - r) * (1 - s)};
    dNdrst[5] = {0.125 * (1 - s) * (1 + t), -0.125 * (1 + r) * (1 + t),
                 0.125 * (1 + r) * (1 - s)};
    dNdrst[6] = {0.125 * (1 + s) * (1 + t), 0.125 * (1 + r) * (1 + t),
                 0.125 * (1 + r) * (1 + s)};
    dNdrst[7] = {-0.125 * (1 + s) * (1 + t), 0.125 * (1 - r) * (1 + t),
                 0.125 * (1 - r) * (1 + s)};
    return dNdrst;
  }
  size_t getNodeCount() const override { return expectedNodeCount; }
};

// 二次ヘキサヘドロン要素
class QuadraticHexahedronElement : public Volume {
 public:
  // static constexpr size_t expectedNodeCount = 20;

  QuadraticHexahedronElement(size_t id, const std::vector<size_t>& nodeIds)
      : Volume(id, GmshElementType::HEXAHEDRON_20NODE, nodeIds, 20) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument(
          "QuadraticHexahedronElement requires exactly 20 nodes");
    }
  }
  size_t getNodeCount() const override { return expectedNodeCount; }
  void computeVolume(std::map<int, Node<double>> nodes) override {
    //   vec3d<double> v1 = nodes[nodeIds[1]].getPos() -
    //   nodes[nodeIds[0]].getPos(); vec3d<double> v2 =
    //   nodes[nodeIds[2]].getPos() - nodes[nodeIds[0]].getPos(); area = 0.5 *
    //   v1.cross(v2).norm();
    throw std::runtime_error(
        "Not implemented yet about compute Volume in "
        "QuadraticHexahedronElement");
  }
  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(20);

    // 8 corner nodes
    N[0] = 0.125 * (1 - r) * (1 - s) * (1 - t) * (-2 - r - s - t);
    N[1] = 0.125 * (1 + r) * (1 - s) * (1 - t) * (-2 + r - s - t);
    N[2] = 0.125 * (1 + r) * (1 + s) * (1 - t) * (-2 + r + s - t);
    N[3] = 0.125 * (1 - r) * (1 + s) * (1 - t) * (-2 - r + s - t);
    N[4] = 0.125 * (1 - r) * (1 - s) * (1 + t) * (-2 - r - s + t);
    N[5] = 0.125 * (1 + r) * (1 - s) * (1 + t) * (-2 + r - s + t);
    N[6] = 0.125 * (1 + r) * (1 + s) * (1 + t) * (-2 + r + s + t);
    N[7] = 0.125 * (1 - r) * (1 + s) * (1 + t) * (-2 - r + s + t);

    // 12 midside nodes
    N[8] = 0.25 * (1 - r * r) * (1 - s) * (1 - t);
    N[9] = 0.25 * (1 + r) * (1 - s * s) * (1 - t);
    N[10] = 0.25 * (1 - r * r) * (1 + s) * (1 - t);
    N[11] = 0.25 * (1 - r) * (1 - s * s) * (1 - t);
    N[12] = 0.25 * (1 - r) * (1 - s) * (1 - t * t);
    N[13] = 0.25 * (1 + r) * (1 - s) * (1 - t * t);
    N[14] = 0.25 * (1 + r) * (1 + s) * (1 - t * t);
    N[15] = 0.25 * (1 - r) * (1 + s) * (1 - t * t);
    N[16] = 0.25 * (1 - r * r) * (1 - s) * (1 + t);
    N[17] = 0.25 * (1 + r) * (1 - s * s) * (1 + t);
    N[18] = 0.25 * (1 - r * r) * (1 + s) * (1 + t);
    N[19] = 0.25 * (1 - r) * (1 - s * s) * (1 + t);

    return N;
  }

  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(20, std::vector<double>(3));

    // 8 corner nodes derivatives
    dNdrst[0] = {-0.125 * (1 - s) * (1 - t) * (-2 - r - s - t) +
                     0.125 * (1 - s) * (1 - t),
                 -0.125 * (1 - r) * (1 - t) * (-2 - r - s - t) +
                     0.125 * (1 - r) * (1 - t),
                 -0.125 * (1 - r) * (1 - s) * (-2 - r - s - t) +
                     0.125 * (1 - r) * (1 - s)};
    dNdrst[1] = {0.125 * (1 - s) * (1 - t) * (-2 + r - s - t) +
                     0.125 * (1 - s) * (1 - t),
                 -0.125 * (1 + r) * (1 - t) * (-2 + r - s - t) +
                     0.125 * (1 + r) * (1 - t),
                 -0.125 * (1 + r) * (1 - s) * (-2 + r - s - t) +
                     0.125 * (1 + r) * (1 - s)};
    dNdrst[2] = {0.125 * (1 + s) * (1 - t) * (-2 + r + s - t) +
                     0.125 * (1 + s) * (1 - t),
                 0.125 * (1 + r) * (1 - t) * (-2 + r + s - t) +
                     0.125 * (1 + r) * (1 - t),
                 -0.125 * (1 + r) * (1 + s) * (-2 + r + s - t) +
                     0.125 * (1 + r) * (1 + s)};
    dNdrst[3] = {-0.125 * (1 + s) * (1 - t) * (-2 - r + s - t) +
                     0.125 * (1 + s) * (1 - t),
                 0.125 * (1 - r) * (1 - t) * (-2 - r + s - t) +
                     0.125 * (1 - r) * (1 - t),
                 -0.125 * (1 - r) * (1 + s) * (-2 - r + s - t) +
                     0.125 * (1 - r) * (1 + s)};
    dNdrst[4] = {-0.125 * (1 - s) * (1 + t) * (-2 - r - s + t) +
                     0.125 * (1 - s) * (1 + t),
                 -0.125 * (1 - r) * (1 + t) * (-2 - r - s + t) +
                     0.125 * (1 - r) * (1 + t),
                 0.125 * (1 - r) * (1 - s) * (-2 - r - s + t) +
                     0.125 * (1 - r) * (1 - s)};
    dNdrst[5] = {0.125 * (1 - s) * (1 + t) * (-2 + r - s + t) +
                     0.125 * (1 - s) * (1 + t),
                 -0.125 * (1 + r) * (1 + t) * (-2 + r - s + t) +
                     0.125 * (1 + r) * (1 + t),
                 0.125 * (1 + r) * (1 - s) * (-2 + r - s + t) +
                     0.125 * (1 + r) * (1 - s)};
    dNdrst[6] = {0.125 * (1 + s) * (1 + t) * (-2 + r + s + t) +
                     0.125 * (1 + s) * (1 + t),
                 0.125 * (1 + r) * (1 + t) * (-2 + r + s + t) +
                     0.125 * (1 + r) * (1 + t),
                 0.125 * (1 + r) * (1 + s) * (-2 + r + s + t) +
                     0.125 * (1 + r) * (1 + s)};
    dNdrst[7] = {-0.125 * (1 + s) * (1 + t) * (-2 - r + s + t) +
                     0.125 * (1 + s) * (1 + t),
                 0.125 * (1 - r) * (1 + t) * (-2 - r + s + t) +
                     0.125 * (1 - r) * (1 + t),
                 0.125 * (1 - r) * (1 + s) * (-2 - r + s + t) +
                     0.125 * (1 - r) * (1 + s)};

    // 12 midside nodes derivatives
    dNdrst[8] = {-0.5 * r * (1 - s) * (1 - t), -0.25 * (1 - r * r) * (1 - t),
                 -0.25 * (1 - r * r) * (1 - s)};
    dNdrst[9] = {0.25 * (1 - s * s) * (1 - t), -0.5 * r * (1 - s * s) * (1 - t),
                 -0.25 * (1 + r) * (1 - s * s)};
    dNdrst[10] = {-0.5 * r * (1 + s) * (1 - t), 0.25 * (1 - r * r) * (1 - t),
                  -0.25 * (1 - r * r) * (1 + s)};
    dNdrst[11] = {-0.25 * (1 - s * s) * (1 - t),
                  -0.5 * r * (1 - s * s) * (1 - t),
                  -0.25 * (1 - r) * (1 - s * s)};
    dNdrst[12] = {-0.25 * (1 - s) * (1 - t * t), -0.25 * (1 - r) * (1 - t * t),
                  -0.5 * t * (1 - r) * (1 - s)};
    dNdrst[13] = {0.25 * (1 - s) * (1 - t * t), -0.25 * (1 + r) * (1 - t * t),
                  -0.5 * t * (1 + r) * (1 - s)};
    dNdrst[14] = {0.25 * (1 + s) * (1 - t * t), 0.25 * (1 + r) * (1 - t * t),
                  -0.5 * t * (1 + r) * (1 + s)};
    dNdrst[15] = {-0.25 * (1 + s) * (1 - t * t), 0.25 * (1 - r) * (1 - t * t),
                  -0.5 * t * (1 - r) * (1 + s)};
    dNdrst[16] = {-0.5 * r * (1 - s) * (1 + t), -0.25 * (1 - r * r) * (1 + t),
                  0.25 * (1 - r * r) * (1 - s)};
    dNdrst[17] = {0.25 * (1 - s * s) * (1 + t),
                  -0.5 * r * (1 - s * s) * (1 + t),
                  0.25 * (1 + r) * (1 - s * s)};
    dNdrst[18] = {-0.5 * r * (1 + s) * (1 + t), 0.25 * (1 - r * r) * (1 + t),
                  0.25 * (1 - r * r) * (1 + s)};
    dNdrst[19] = {-0.25 * (1 - s * s) * (1 + t),
                  -0.5 * r * (1 - s * s) * (1 + t),
                  0.25 * (1 - r) * (1 - s * s)};

    return dNdrst;
  }
};

// // ピラミッド要素
class PyramidElement : public Volume {
 public:
  // static constexpr size_t expectedNodeCount = 5;

  PyramidElement(size_t id, const std::vector<size_t>& nodeIds)
      : Volume(id, GmshElementType::PYRAMID_5NODE, nodeIds, 5) {
    if (nodeIds.size() != expectedNodeCount) {
      throw std::invalid_argument("PyramidElement requires exactly 5 nodes");
    }
  }

  size_t getNodeCount() const override { return expectedNodeCount; }
  void computeVolume(std::map<int, Node<double>> nodes) override {
    //   vec3d<double> v1 = nodes[nodeIds[1]].getPos() -
    //   nodes[nodeIds[0]].getPos(); vec3d<double> v2 =
    //   nodes[nodeIds[2]].getPos() - nodes[nodeIds[0]].getPos(); area = 0.5 *
    //   v1.cross(v2).norm();
    throw std::runtime_error(
        "Not implemented yet about compute Volume in "
        "PyramidElement");
  }

  std::vector<double> computeShapeFunctions(const double& r, const double& s,
                                            const double& t) const override {
    std::vector<double> N(5);

    N[0] = (1 - r - s - t) * (1 + r - s - t) * (1 - r + s - t) *
           (1 + r + s - t) * (1 - t) / 16;
    N[1] = (1 + r - s - t) * (1 - r + s - t) * (1 + r + s - t) * (1 - t) / 8;
    N[2] = (1 - r + s - t) * (1 + r + s - t) * (1 - t) / 8;
    N[3] = (1 + r + s - t) * (1 - t) / 4;
    N[4] = t;

    return N;
  }
  std::vector<std::vector<double>> computeShapeFunctionsDerivative(
      const double& r, const double& s, const double& t) const override {
    std::vector<std::vector<double>> dNdrst(5, std::vector<double>(3));

    // 形状関数の導関数の計算
    dNdrst[0] = {
        -0.25 * (1 - s - t) * (1 - 2 * r + s + t) * (1 - t),
        -0.25 * (1 - r - t) * (1 + r - 2 * s + t) * (1 - t),
        -0.25 * (1 - r - s) * (1 + r - s) * (1 - r + s) * (1 + r + s) / 16 -
            0.25 * (1 - r - s - t) * (1 + r - s - t) * (1 - r + s - t) *
                (1 + r + s - t)};

    dNdrst[1] = {0.25 * (1 - s - t) * (1 + r - s - t) * (1 - t),
                 -0.25 * (1 + r - t) * (1 - 2 * s + r + t) * (1 - t),
                 -0.25 * (1 + r - s) * (1 - s) * (1 + s) * (1 - t) / 8};

    dNdrst[2] = {-0.25 * (1 + s - t) * (1 - 2 * r + s + t) * (1 - t),
                 0.25 * (1 - r - t) * (1 + r + s - t) * (1 - t),
                 -0.25 * (1 - r + s) * (1 - r) * (1 + r) * (1 - t) / 8};

    dNdrst[3] = {
        0.25 * (1 + s - t) * (1 + r + s - t) * (1 - t),
        0.25 * (1 + r - t) * (1 + r + s - t) * (1 - t),
        -0.25 * (1 + r + s) * (1 - r) * (1 + r) * (1 - s) * (1 + s) / 4};

    dNdrst[4] = {0, 0, 1};

    return dNdrst;
  }
};

// // 二次ピラミッド要素
// class QuadraticPyramidElement : public Volume {
//  public:
//   static constexpr size_t expectedNodeCount = 13;

//   QuadraticPyramidElement(size_t id, const std::vector<size_t>& nodeIds)
//       : Volume(id, GmshElementType::PYRAMID_13NODE, nodeIds) {
//     if (nodeIds.size() != expectedNodeCount) {
//       throw std::invalid_argument(
//           "QuadraticPyramidElement requires exactly 13 nodes");
//     }
//   }

//   size_t getNodeCount() const override { return expectedNodeCount; }
// };

class Mesh {
 public:
  Mesh() = default;
  Mesh(Mesh&& other) noexcept
      : nodes(std::move(other.nodes)),
        elements(std::move(other.elements)),
        physicalGrMap(std::move(other.physicalGrMap)),
        physicalGrToElements(std::move(other.physicalGrToElements)) {}
  Mesh& operator=(Mesh&& other) noexcept {
    if (this != &other) {
      nodes = std::move(other.nodes);
      elements = std::move(other.elements);
      physicalGrMap = std::move(other.physicalGrMap);
      physicalGrToElements = std::move(other.physicalGrToElements);
    }
    return *this;
  }
  std::map<int, Node<double>> nodes;
  std::vector<std::shared_ptr<Element>> elements;
  // surface名とsurfaceのIDをマッピング
  std::unordered_map<std::string, int> physicalGrMap;
  // IDをキーとし、そのphysicalGrに所属するelementのIDリストを値とするマップ
  std::unordered_map<int, std::vector<std::shared_ptr<Element>>>
      physicalGrToElements;

  /// volumeのIDとそのvolumeに対応するmaterialのIDをペアで格納
  std::vector<std::pair<int, int>> volumeForMaterial;

  ///
  void addNode(const Node<double>& node) { nodes[node.getId()] = node; }

  void addElement(std::shared_ptr<Element>& element) {
    elements.push_back(std::move(element));
  }
  const Node<double>& getNode(size_t nodeId) const { return nodes.at(nodeId); }

  // const std::unique_ptr<Element>& getElement(size_t elementIndex) const {
  const std::shared_ptr<Element>& getElement(size_t elementIndex) const {
    if (elementIndex < elements.size()) {
      return elements[elementIndex];
    } else {
      throw std::out_of_range("Element index is out of range");
    }
  }

  size_t getNumberOfNodes() const { return nodes.size(); }
  size_t getNumberOfElements() const { return elements.size(); }
  size_t getNumberOfSurfaces(const CFDConfig& config) const {
    size_t sum(0);
    for (auto& it : config.boundaryConditions) {
      sum += physicalGrToElements.at(it.id).size();
    }
    if (physicalGrMap.find("default") != physicalGrMap.end())
      sum += physicalGrToElements.at(physicalGrMap.at("default")).size();
    return sum;
  }
  size_t getNumberOfVolumes() const {
    size_t sum(0);
    for (auto& it : volumeForMaterial) {
      sum += physicalGrToElements.at(it.first).size();
    }
    return sum;
  }

  /// @ brief set Pair of Physcial Gropuup ID and Material ID
  void setVolumeForMaterial(std::string volumename, const int id) {
    for (auto& it : physicalGrMap) {
      if (it.first == volumename) {
        volumeForMaterial.push_back(std::make_pair(it.second, id));
      }
    }
  }
};

}  // namespace FemFy
#endif  // _MESH_H_
