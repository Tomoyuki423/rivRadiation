// Created by Tomoyuki Kawashima on 27/June/2023
#ifndef _READMSH_H_
#define _READMSH_H_

#include <memory>
#ifdef _WIN32
#include <gmsh-sdk-win/include/gmsh.h>
#elif defined(__linux__)
#include <gmsh-sdk-linux/include/gmsh.h>
#else
#error "Unsupported operating system"
#endif

#include <mesh.h>

#include <iostream>
#include <typeinfo>
#include <vector>

namespace FemFy {
// template <typename T>
// T checkTypeMatch(int type, size_t id, const std::vector<size_t>& nodeIds) {
//   switch (type) {
//     case static_cast<int>(GmshElementType::LINE_2NODE): {
//       LineElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::LINE_3NODE): {
//       QuadraticLineElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::TRIANGLE_3NODE): {
//       TriangleElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::TETRAHEDRON_4NODE): {
//       TetrahedronElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::QUAD_4NODE): {
//       QuadrangleElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::TRIANGLE_6NODE): {
//       QuadraticTriangleElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::QUAD_8NODE): {
//       QuadraticQuadrangleElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::PRISM_6NODE): {
//       PrismElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::PRISM_15NODE): {
//       QuadraticPrismElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::HEXAHEDRON_8NODE): {
//       HexahedronElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::HEXAHEDRON_20NODE): {
//       QuadraticHexahedronElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::PYRAMID_5NODE): {
//       PyramidElement element(id, nodeIds);
//       return element;
//     }
//     case static_cast<int>(GmshElementType::PYRAMID_13NODE): {
//       QuadraticPrismElement element(id, nodeIds);
//       return element;
//     }
//     default:
//       return nullptr;
//   }
// }
std::unique_ptr<Element> checkTypeMatch(int type, size_t id,
                                        const std::vector<size_t>& nodeIds) {
  switch (type) {
    case static_cast<int>(GmshElementType::POINT_1NODE):
      return std::make_unique<Point>(id, nodeIds);
    case static_cast<int>(GmshElementType::LINE_2NODE):
      return std::make_unique<LineElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::LINE_3NODE):
      return std::make_unique<QuadraticLineElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::TRIANGLE_3NODE):
      return std::make_unique<TriangleElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::TETRAHEDRON_4NODE):
      return std::make_unique<TetrahedronElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::QUAD_4NODE):
      return std::make_unique<QuadrangleElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::TRIANGLE_6NODE):
      return std::make_unique<QuadraticTriangleElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::QUAD_8NODE):
      return std::make_unique<QuadraticQuadrangleElement>(id, nodeIds);
    // case static_cast<int>(GmshElementType::PRISM_6NODE):
    //   return std::make_unique<PrismElement>(id, nodeIds);
    // case static_cast<int>(GmshElementType::PRISM_15NODE):
    //   return std::make_unique<QuadraticPrismElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::HEXAHEDRON_8NODE):
      return std::make_unique<HexahedronElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::HEXAHEDRON_20NODE):
      return std::make_unique<QuadraticHexahedronElement>(id, nodeIds);
    // case static_cast<int>(GmshElementType::PYRAMID_5NODE):
    //   return std::make_unique<PyramidElement>(id, nodeIds);
    // case static_cast<int>(GmshElementType::PYRAMID_13NODE):
    //   return std::make_unique<QuadraticPyramidElement>(id, nodeIds);
    default:
      return nullptr;
  }
}

void readMSHFile(const std::string& fileName, Mesh& mesh) {
  gmsh::initialize();
  gmsh::open(fileName);

  std::vector<size_t> nodeTags;
  std::vector<double> coord, parametricCoord;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

  for (size_t i = 0; i < nodeTags.size(); i++) {
    Node<double> node(nodeTags[i], coord[3 * i], coord[3 * i + 1],
                      coord[3 * i + 2]);
    mesh.addNode(node);
  }

  std::vector<int> elementTypes;
  gmsh::model::mesh::getElementTypes(elementTypes);

  for (auto type : elementTypes) {
    std::vector<size_t> elementTags, nodeTags;
    gmsh::model::mesh::getElementsByType(type, elementTags, nodeTags);
    size_t numNodes = nodeTags.size() / elementTags.size();
    std::vector<size_t> nodeTagsTempo(numNodes);
    for (size_t i = 0; i < elementTags.size(); i++) {
      for (size_t j = 0; j < numNodes; j++) {
        nodeTagsTempo[j] = nodeTags[i * numNodes + j];
      }
      std::unique_ptr<Element> element =
          checkTypeMatch(type, elementTags[i], nodeTagsTempo);
      mesh.addElement(element);
      //     elements.push_back(element);
    }
  }

  gmsh::finalize();
}
}  // namespace FemFy
#endif  // _READMSH_H_
