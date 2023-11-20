#include <CFDConfig.h>
// #include <Simulation.h>
#include <mesh.h>
#include <readmsh.h>

#include <iostream>
#include <string>

using namespace std;
using namespace FemFy;

int main() {
  cout << "-----------Program start--------------" << endl;
  string configFileName = "./config.json";
  string meshFileName = "./msh/cube-Body.msh";
  CFDConfig config;
  Mesh mesh;
  config.loadFromConfigFile(configFileName);
  readMSHFile(meshFileName, mesh);
  int i(0);
  std::cout << "mesh.elements.size() = " << mesh.elements.size() << std::endl;
  for (auto& it : mesh.elements) {
    cout << i++ << endl;
    it->printInfo();
  }
  //   Simulation sim(config, mesh);
  //   sim.run();
  cout << "Finish sim" << endl;
  return 0;
}
