// Created by Tomoyuki Kawashima on 27/June/2023
// include guard
#ifndef _READMSH_H_
#define _READMSH_H_

#include <gmsh.h>

#include <iostream>
#include <vector>

using namespace std;
namespace IM {

void readMSHFile(const string& fileName, vector<Node>& nodes,
                 vector<Element>& elements) {
  gmsh::initialize();
  gmsh::open(fileName);

  vector<size_t> nodeTags;
  vector<double> coord, parametricCoord;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

  for (size_t i = 0; i < nodeTags.size(); i++) {
    Node node;
    node.id = nodeTags[i];
    node.x = coord[3 * i];
    node.y = coord[3 * i + 1];
    node.z = coord[3 * i + 2];
    nodes.push_back(node);
  }

  vector<int> elementTypes;
  gmsh::model::mesh::getElementTypes(elementTypes);

  for (auto type : elementTypes) {
    vector<size_t> elementTags, nodeTags;
    gmsh::model::mesh::getElementsByType(type, elementTags, nodeTags);
    size_t numNodes = nodeTags.size() / elementTags.size();
    for (size_t i = 0; i < elementTags.size(); i++) {
      Element element;
      element.id = elementTags[i];
      element.type = type;
      for (size_t j = 0; j < numNodes; j++) {
        element.nodeid.push_back(nodeTags[i * numNodes + j]);
      }
      elements.push_back(element);
    }
  }

  gmsh::finalize();
}
}  // namespace IM
#endif  // _READMSH_H_
