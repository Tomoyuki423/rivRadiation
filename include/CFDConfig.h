#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <fstream>
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <vector>

using json = nlohmann::json;
namespace FemFy {
struct GridSize {
  int x, y, z;
};

struct BoundaryCondition {
  std::string name;
  std::string type;
  std::string face;
  std::map<std::string, double> value;
};

class CFDConfig {
 public:
  std::string solver;
  GridSize gridSize;
  double timeStep;
  int iterations;
  std::vector<BoundaryCondition> boundaryConditions;
  bool isTimeDependency;

  void loadFromConfigFile(const std::string& filepath) {
    std::ifstream inputFile(filepath);
    if (!inputFile) {
      std::cerr << "Could not open file: " << filepath << std::endl;
      return;
    }

    json j;
    inputFile >> j;

    if (j.contains("solver")) solver = j["solver"];
    if (j.contains("timeDependency")) {
      std::string timeDependency = j["timeDependency"];
      if (timeDependency == "steady") {
        // 定常解析の設定
        isTimeDependency = false;
      } else if (timeDependency == "transient") {
        // 非定常解析の設定
        isTimeDependency = true;
      }
    }
    if (j.contains("gridSize")) {
      if (j["gridSize"].contains("x")) gridSize.x = j["gridSize"]["x"];
      if (j["gridSize"].contains("y")) gridSize.y = j["gridSize"]["y"];
      if (j["gridSize"].contains("z")) gridSize.z = j["gridSize"]["z"];
    }
    if (j.contains("timeStep")) timeStep = j["timeStep"];
    if (j.contains("iterations")) iterations = j["iterations"];

    if (j.contains("boundaryConditions")) {
      for (const auto& bc : j["boundaryConditions"]) {
        BoundaryCondition boundaryCondition;
        if (bc.contains("name")) boundaryCondition.name = bc["name"];
        if (bc.contains("type")) boundaryCondition.type = bc["type"];
        if (bc.contains("face")) boundaryCondition.face = bc["face"];
        if (bc.contains("value"))
          boundaryCondition.value =
              bc["value"].get<std::map<std::string, double>>();
        boundaryConditions.push_back(boundaryCondition);
      }
    }
  }
};
}  // namespace FemFy
#endif  // _CONFIG_H_
