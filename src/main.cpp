#include <FemFy/FemFy.hpp>

using namespace std;
using namespace FemFy;

int main() {
  cout << "-----------Program start--------------" << endl;
  string configFileName = "./config.json";
  string meshFileName = "./msh/cube-Body_10_30.msh";
  CFDConfig config;
  Mesh mesh;
  config.loadFromConfigFile(configFileName);
  readMSHFile(meshFileName, mesh);
  Simulation sim(config, std::move(mesh));
  sim.run();
  cout << "Finish sim" << endl;
  return 0;
}
