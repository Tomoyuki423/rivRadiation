// Created by Tomoyuki Kawashima on 27/June/2023
#ifndef _READMSH_H_
#define _READMSH_H_

#include <FemFy/FemFy.hpp>
#include <FemFy/mesh.hpp>
#include <FemFy/utils.hpp>
////////////////////////////////////////////
#include <memory>
#ifdef _WIN32
#include <gmsh-sdk-win/include/gmsh.h>
#elif defined(__linux__)
#include <gmsh-sdk-linux/include/gmsh.h>
#else
#error "Unsupported operating system"
#endif

#include <algorithm>
#include <iostream>
#include <typeinfo>
#include <unordered_map>
#include <vector>
namespace FemFy {
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
    case static_cast<int>(GmshElementType::TETRAHEDRON_10NODE):
      return std::make_unique<QuadraticTetrahedronElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::QUAD_4NODE):
      return std::make_unique<QuadrangleElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::TRIANGLE_6NODE):
      return std::make_unique<QuadraticTriangleElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::QUAD_8NODE):
      return std::make_unique<QuadraticQuadrangleElement>(id, nodeIds);
    case static_cast<int>(GmshElementType::PRISM_6NODE):
      return std::make_unique<PrismElement>(id, nodeIds);
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

bool checkPoint(const Mesh& mesh, const int index) {
  if (mesh.elements[index]->getType() == GmshElementType::POINT_1NODE) {
    return true;
  }
  return false;
}
bool checkLine(const Mesh& mesh, const int index) {
  if (mesh.elements[index]->getType() == GmshElementType::LINE_2NODE ||
      mesh.elements[index]->getType() == GmshElementType::LINE_3NODE) {
    return true;
  }
  return false;
}
bool checkSurface(const Mesh& mesh, const int index) {
  if (mesh.elements[index]->getType() == GmshElementType::TRIANGLE_3NODE ||
      mesh.elements[index]->getType() == GmshElementType::TRIANGLE_6NODE ||
      mesh.elements[index]->getType() == GmshElementType::QUAD_4NODE ||
      mesh.elements[index]->getType() == GmshElementType::QUAD_8NODE ||
      mesh.elements[index]->getType() == GmshElementType::QUAD_9NODE) {
    return true;
  }
  return false;
}
bool checkVolume(const Mesh& mesh, const int index) {
  if (mesh.elements[index]->getType() == GmshElementType::TETRAHEDRON_4NODE ||
      mesh.elements[index]->getType() == GmshElementType::TETRAHEDRON_10NODE ||
      mesh.elements[index]->getType() == GmshElementType::HEXAHEDRON_8NODE ||
      mesh.elements[index]->getType() == GmshElementType::HEXAHEDRON_20NODE ||
      mesh.elements[index]->getType() == GmshElementType::PRISM_6NODE ||
      mesh.elements[index]->getType() == GmshElementType::PRISM_15NODE ||
      mesh.elements[index]->getType() == GmshElementType::PYRAMID_5NODE ||
      mesh.elements[index]->getType() == GmshElementType::PYRAMID_13NODE) {
    return true;
  }
  return false;
}

void readMSHFile(const std::string& fileName, Mesh& mesh) {
  gmsh::initialize();
  gmsh::open(fileName);
  std::string name;
  gmsh::model::getCurrent(name);
  std::cout << "Model " << name << " (" << gmsh::model::getDimension()
            << "D)\n";
  std::vector<std::pair<int, int>> entities;
  gmsh::model::getEntities(entities);
  int countElem(0);
  int countPoint(0);
  int startPoint(0);
  int countLine(0);
  int startLine(0);
  int countSurface(0);
  int startSurface(0);
  int countVolume(0);
  int startVolume(0);
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
      std::shared_ptr<Element> element =
          checkTypeMatch(type, elementTags[i], nodeTagsTempo);
      mesh.addElement(element);
      //     elements.push_back(element);
    }
  }

  for (auto e : entities) {
    int dim = e.first, tag = e.second;
    std::vector<double> nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);
    std::string type;
    gmsh::model::getType(dim, tag, type);
    std::string name;
    gmsh::model::getEntityName(dim, tag, name);
    if (name.size()) name += " ";
    std::cout << "Entity " << name << "(" << dim << "," << tag << ") of type "
              << type << "\n";
    int numElem = 0;
    for (auto& tags : elemTags) numElem += tags.size();
    std::cout << " - Mesh has " << nodeTags.size() << " nodes and " << numElem
              << " elements\n";
    // for (int i = countElem; i < countElem + numElem; i++) { //   std::cout <<
    // "mesh.elements[" << i << "]->printInfo() = ";
    //   mesh.elements[i]->printInfo();
    // }

    // * Does the entity belong to physical groups?
    std::vector<int> physicalTags;
    gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);
    ///
    if (dim == 2 && physicalTags.size() == 0) {
      mesh.physicalGrMap["default"] = 0;  // add by kawashima
      if (type == "Discrete surface") {
        if (countSurface == 0) {
          for (int i = 0, size = mesh.elements.size(); i < size; i++) {
            if (checkSurface(mesh, i)) {
              startSurface = i;
              break;
            }
          }
        }
        for (int i = startSurface + countSurface;
             i < startSurface + countSurface + numElem; i++) {
          // mesh.elements[i]->printInfo();
          mesh.physicalGrToElements[mesh.physicalGrMap["default"]].push_back(
              mesh.elements[i]);
        }
        countSurface += numElem;
      }
    }
    ///
    if (dim == 3 && physicalTags.size() == 0) {
      std::cerr << "Physical group for volume don't be not defined"
                << std::endl;
      exit(1);
    }
    if (physicalTags.size()) {
      std::cout << " - Physical group: ";
      for (auto physTag : physicalTags) {
        std::string n;
        std::string nn;
        gmsh::model::getPhysicalName(dim, physTag, n);
        mesh.physicalGrMap[n] = physTag;  // add by kawashima
        if (n.size()) nn = n + " ";
        std::cout << nn << "(" << dim << ", " << physTag << ") ";

        if (type == "Discrete point") {
          if (countPoint == 0) {
            for (int i = 0, size = mesh.elements.size(); i < size; i++) {
              if (checkPoint(mesh, i)) {
                startPoint = i;
                break;
              }
            }
          }
          for (int i = startPoint + countPoint;
               i < startPoint + countPoint + numElem; i++) {
            // mesh.elements[i]->printInfo();
            mesh.physicalGrToElements[mesh.physicalGrMap[n]].push_back(
                mesh.elements[i]);
          }
          countPoint += numElem;
        } else if (type == "Discrete curve") {
          if (countSurface == 0) {
            for (int i = 0, size = mesh.elements.size(); i < size; i++) {
              if (checkLine(mesh, i)) {
                startLine = i;
                break;
              }
            }
          }
          for (int i = startLine + countLine;
               i < startLine + countLine + numElem; i++) {
            // mesh.elements[i]->printInfo();
            mesh.physicalGrToElements[mesh.physicalGrMap[n]].push_back(
                mesh.elements[i]);
          }
          countLine += numElem;

        } else if (type == "Discrete surface") {
          if (countSurface == 0) {
            for (int i = 0, size = mesh.elements.size(); i < size; i++) {
              if (checkSurface(mesh, i)) {
                startSurface = i;
                break;
              }
            }
          }
          for (int i = startSurface + countSurface;
               i < startSurface + countSurface + numElem; i++) {
            // mesh.elements[i]->printInfo();
            mesh.physicalGrToElements[mesh.physicalGrMap[n]].push_back(
                mesh.elements[i]);
          }
          countSurface += numElem;
        } else if (type == "Discrete volume") {
          if (countVolume == 0) {
            for (int i = 0, size = mesh.elements.size(); i < size; i++) {
              if (checkVolume(mesh, i)) {
                startVolume = i;
                break;
              }
            }
          }
          for (int i = startVolume + countVolume;
               i < startVolume + countVolume + numElem; i++) {
            // mesh.elements[i]->printInfo();
            mesh.physicalGrToElements[mesh.physicalGrMap[n]].push_back(
                mesh.elements[i]);
          }
          countVolume += numElem;
        } else {
          std::cout
              << "Unknown type of entity. Please check the type of entity."
              << type << std::endl;
        }
      }
      std::cout << "\n";
    }
    // physicalGrToElements[physicalGrMap[n]].push_back(
    //     mesh.elements[countElem + i]);
    countElem += numElem;
  }

  // We can use this to clear all the model data:
  gmsh::clear();

  gmsh::finalize();
}

// void test2(const std::string& fileName, Mesh& mesh) {
//   gmsh::initialize();
//   gmsh::open(fileName);
//   std::string name;
//   gmsh::model::getCurrent(name);
//   std::cout << "Model " << name << " (" << gmsh::model::getDimension()
//             << "D)\n";
//   std::vector<std::pair<int, int>> entities;
//   gmsh::model::getEntities(entities);
//   int countElem(0);
//   int countPoint(0);
//   int countLine(0);
//   int startLine(0);
//   int countSurface(0);
//   int startSurface(0);
//   int countVolume(0);
//   int startVolume(0);

//   ////////////////////////////////////
//   // surface名とsurfaceのIDをマッピング
//   std::unordered_map<std::string, int> physicalGrMap;
//   // IDをキーとし、そのsurfaceに所属するelementのIDリストを値とするマップ
//   std::unordered_map<int, std::vector<std::shared_ptr<Element>>>
//       physicalGrToElements;
//   ////////////////////////////////////
//   for (auto e : entities) {
//     int dim = e.first, tag = e.second;
//     // if (dim == 0 || dim == 1) continue;
//     std::vector<std::size_t> nodeTags;
//     std::vector<double> nodeCoords, nodeParams;
//     gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim,
//     tag); std::vector<int> elemTypes; std::vector<std::vector<std::size_t>>
//     elemTags, elemNodeTags; gmsh::model::mesh::getElements(elemTypes,
//     elemTags, elemNodeTags, dim, tag); std::string type;
//     gmsh::model::getType(dim, tag, type); std::string name;
//     gmsh::model::getEntityName(dim, tag, name); if (name.size()) name
//     += " "; std::cout << "Entity " << name << "(" << dim << "," << tag <<
//     ") of type "
//               << type << "\n";
//     int numElem = 0;
//     std::cout << "elemTags.size() = " << elemTags.size() << std::endl;
//     for (auto& tags : elemTags) numElem += tags.size();
//     std::cout << " - Mesh has " << nodeTags.size() << " nodes and " <<
//     numElem
//               << " elements\n";
//     // for (int i = countElem; i < countElem + numElem; i++) {
//     //   std::cout << "mesh.elements[" << i << "]->printInfo() = ";
//     //   mesh.elements[i]->printInfo();
//     // }

//     // * Does the entity belong to physical groups?
//     std::vector<int> physicalTags;
//     gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);
//     std::cout << "physicalTags.size() = " << physicalTags.size() <<
//     std::endl; if (physicalTags.size()) {
//       std::cout << " - Physical group: ";
//       for (auto physTag : physicalTags) {
//         std::string n;
//         std::string nn;
//         gmsh::model::getPhysicalName(dim, physTag, n);
//         physicalGrMap[n] = physTag;  // add by kawashima
//         if (n.size()) nn = n + " ";
//         std::cout << nn << "(" << dim << ", " << physTag << ") ";

//         if (type == "Discrete point") {
//         } else if (type == "Discrete curve") {
//           if (countSurface == 0) {
//             for (int i = 0, size = mesh.elements.size(); i < size; i++) {
//               if (checkLine(mesh, i)) {
//                 startLine = i;
//                 break;
//               }
//             }
//           }
//           for (int i = startLine + countLine;
//                i < startLine + countLine + numElem; i++) {
//             std::cout << "mesh.elements[" << i << "]->printInfo() = ";
//             mesh.elements[i]->printInfo();
//             physicalGrToElements[physicalGrMap[n]].push_back(mesh.elements[i]);
//           }
//           countLine += numElem;

//         } else if (type == "Discrete surface") {
//           if (countSurface == 0) {
//             for (int i = 0, size = mesh.elements.size(); i < size; i++) {
//               if (checkSurface(mesh, i)) {
//                 startSurface = i;
//                 break;
//               }
//             }
//           }
//           for (int i = startSurface + countSurface;
//                i < startSurface + countSurface + numElem; i++) {
//             std::cout << "mesh.elements[" << i << "]->printInfo() = ";
//             mesh.elements[i]->printInfo();
//             physicalGrToElements[physicalGrMap[n]].push_back(mesh.elements[i]);
//           }
//           countSurface += numElem;
//         } else if (type == "Discrete volume") {
//           if (countVolume == 0) {
//             for (int i = 0, size = mesh.elements.size(); i < size; i++) {
//               if (checkVolume(mesh, i)) {
//                 startVolume = i;
//                 break;
//               }
//             }
//           }
//           for (int i = startVolume + countVolume;
//                i < startVolume + countVolume + numElem; i++) {
//             std::cout << "mesh.elements[" << i << "]->printInfo() = ";
//             mesh.elements[i]->printInfo();
//             physicalGrToElements[physicalGrMap[n]].push_back(mesh.elements[i]);
//           }
//           countVolume += numElem;
//         } else {
//           std::cout
//               << "Unknown type of entity. Please check the type of entity."
//               << type << std::endl;
//         }
//       }
//       std::cout << "\n";
//     }
//     // physicalGrToElements[physicalGrMap[n]].push_back(
//     //     mesh.elements[countElem + i]);
//     countElem += numElem;
//   }
//   for (auto& it : physicalGrMap) {
//     std::cout << "physicalGrMap[" << it.first << "] = " << it.second
//               << std::endl;
//   }
//   for (auto& it : physicalGrToElements[physicalGrMap["volume1"]]) {
//     it->printInfo();
//   }
//   // We can use this to clear all the model data:
//   gmsh::clear();

//   gmsh::finalize();
//   return;
// }

// void test(const std::string& fileName) {
//   gmsh::initialize();

//   gmsh::open(fileName);

//   // Print the model name and dimension:
//   std::string name;
//   gmsh::model::getCurrent(name);
//   std::cout << "Model " << name << " (" << gmsh::model::getDimension()
//             << "D)\n";

//   // Geometrical data is made of elementary model `entities', called
//   `points'
//   // (entities of dimension 0), `curves' (entities of dimension 1),
//   `surfaces'
//   // (entities of dimension 2) and `volumes' (entities of dimension 3). As
//   we
//   // have seen in the other C++ tutorials, elementary model entities are
//   // identified by their dimension and by a `tag': a strictly positive
//   // identification number. Model entities can be either CAD entities (from
//   the
//   // built-in `geo' kernel or from the OpenCASCADE `occ' kernel) or
//   `discrete'
//   // entities (defined by a mesh). `Physical groups' are collections of
//   model
//   // entities and are also identified by their dimension and by a tag.

//   // Get all the elementary entities in the model, as a vector of
//   (dimension,
//   // tag) pairs:
//   std::vector<std::pair<int, int>> entities;
//   gmsh::model::getEntities(entities);

//   for (auto e : entities) {
//     // Dimension and tag of the entity:
//     int dim = e.first, tag = e.second;

//     // Mesh data is made of `elements' (points, lines, triangles, ...),
//     defined
//     // by an ordered list of their `nodes'. Elements and nodes are
//     identified by
//     // `tags' as well (strictly positive identification numbers), and are
//     stored
//     // ("classified") in the model entity they discretize. Tags for
//     elements and
//     // nodes are globally unique (and not only per dimension, like
//     entities).

//     // A model entity of dimension 0 (a geometrical point) will contain a
//     mesh
//     // element of type point, as well as a mesh node. A model curve will
//     contain
//     // line elements as well as its interior nodes, while its boundary
//     nodes
//     // will be stored in the bounding model points. A model surface will
//     contain
//     // triangular and/or quadrangular elements and all the nodes not
//     classified
//     // on its boundary or on its embedded entities. A model volume will
//     contain
//     // tetrahedra, hexahedra, etc. and all the nodes not classified on its
//     // boundary or on its embedded entities.

//     // Get the mesh nodes for the entity (dim, tag):
//     std::vector<std::size_t> nodeTags;
//     std::vector<double> nodeCoords, nodeParams;
//     gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim,
//     tag);

//     // Get the mesh elements for the entity (dim, tag):
//     std::vector<int> elemTypes;
//     std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
//     gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim,
//     tag);

//     // Elements can also be obtained by type, by using `getElementTypes()'
//     // followed by `getElementsByType()'.

//     // Let's print a summary of the information available on the entity and
//     its
//     // mesh.

//     // * Type of the entity:
//     std::string type;
//     gmsh::model::getType(dim, tag, type);
//     std::string name;
//     gmsh::model::getEntityName(dim, tag, name);
//     if (name.size()) name += " ";
//     std::cout << "Entity " << name << "(" << dim << "," << tag << ") of
//     type
//     "
//               << type << "\n";

//     // * Number of mesh nodes and elements:
//     int numElem = 0;
//     for (auto& tags : elemTags) numElem += tags.size();
//     std::cout << " - Mesh has " << nodeTags.size() << " nodes and " <<
//     numElem
//               << " elements\n";

//     // * Upward and downward adjacencies:
//     std::vector<int> up, down;
//     gmsh::model::getAdjacencies(dim, tag, up, down);
//     if (up.size()) {
//       std::cout << " - Upward adjacencies: ";
//       for (auto e : up) std::cout << e << " ";
//       std::cout << "\n";
//     }
//     if (down.size()) {
//       std::cout << " - Downward adjacencies: ";
//       for (auto e : down) std::cout << e << " ";
//       std::cout << "\n";
//     }

//     // * Does the entity belong to physical groups?
//     std::vector<int> physicalTags;
//     gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);
//     if (physicalTags.size()) {
//       std::cout << " - Physical group: ";
//       for (auto physTag : physicalTags) {
//         std::string n;
//         gmsh::model::getPhysicalName(dim, physTag, n);
//         if (n.size()) n += " ";
//         std::cout << n << "(" << dim << ", " << physTag << ") ";
//       }
//       std::cout << "\n";
//     }

//     // * Is the entity a partition entity? If so, what is its parent
//     entity? std::vector<int> partitions; gmsh::model::getPartitions(dim,
//     tag, partitions); if (partitions.size()) {
//       std::cout << " - Partition tags:";
//       for (auto part : partitions) std::cout << " " << part;
//       int parentDim, parentTag;
//       gmsh::model::getParent(dim, tag, parentDim, parentTag);
//       std::cout << " - parent entity (" << parentDim << "," << parentTag
//                 << ")\n";
//     }

//     // * List all types of elements making up the mesh of the entity:
//     for (auto elemType : elemTypes) {
//       std::string name;
//       int d, order, numv, numpv;
//       std::vector<double> param;
//       gmsh::model::mesh::getElementProperties(elemType, name, d, order,
//       numv,
//                                               param, numpv);
//       std::cout << " - Element type: " << name << ", order " << order <<
//       "\n"; std::cout << "   with " << numv << " nodes in param coord: (";
//       for (auto p : param) std::cout << p << " ";
//       std::cout << ")\n";
//     }
//   }

//   // We can use this to clear all the model data:
//   gmsh::clear();

//   gmsh::finalize();
//   return;
// }

}  // namespace FemFy
#endif  // _READMSH_H_
