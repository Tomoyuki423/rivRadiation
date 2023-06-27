// Compile: g++ -std=c++11 -o readMSH readMSH.cpp
// Run: ./readMSH
// Output: 0.5

#include <element.h>
#include <readmsh.h>
#include <vector3d.h>

#include <iostream>
#include <vector>

int main() {
  vector<Node> nodes;
  vector<Element> elements;
  readMSHFile("./msh/mesh.msh", nodes, elements);
  // Now nodes and elements are filled with data from the .msh file
  cout << nodes.size() << endl;
  for (auto node : nodes) {
    cout << node.x << endl;
  }
  return 0;
}
