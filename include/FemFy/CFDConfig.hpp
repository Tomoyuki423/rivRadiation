#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <FemFy/FemFy.hpp>
////////////////////////////////
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <stdexcept>

using json = nlohmann::json;
namespace FemFy {
struct GridSize {
  int NbNode;
  int NbLine;
  int NbSurface;
  int NbVolume;
};

struct BoundaryCondition {
  int id;
  std::string name;
  std::string kind;
  std::string type;
  double value;
  double flux;
};

struct VolumeProperty {
  int id;
  std::string name;
  std::string kind;
  // std::string type;
  std::string material;
};

struct MaterialProperty {
  int id;
  std::string name;
  std::string kind;
  std::string type;
  double density;
  double specificHeat;
  double thermalConductivity;
};

class CFDConfig {
 public:
  std::string solver;
  GridSize gridSize;
  double timeStep;
  int iterations;
  std::vector<BoundaryCondition> boundaryConditions;
  std::vector<VolumeProperty> volumeProperties;
  std::vector<MaterialProperty> materials;
  bool isTimeDependency;
  inline void loadFromConfigFile(const std::string& filepath);
};

void CFDConfig::loadFromConfigFile(const std::string& filepath) {
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
  // if (j.contains("gridSize")) {
  //   if (j["gridSize"].contains("x")) gridSize.x = j["gridSize"]["x"];
  //   if (j["gridSize"].contains("y")) gridSize.y = j["gridSize"]["y"];
  //   if (j["gridSize"].contains("z")) gridSize.z = j["gridSize"]["z"];
  // }
  if (j.contains("timeStep")) timeStep = j["timeStep"];
  if (j.contains("iterations")) iterations = j["iterations"];

  bool isDirichlet(false);
  bool isNuemann(false);
  bool isThird(false);
  bool isName(false);
  bool isKind(false);
  bool isType(false);
  bool isFlux(false);
  bool isValue(false);
  bool isMaterial(false);
  bool isDensity(false);
  bool isSpecificHeat(false);
  bool isThermalConductivity(false);
  bool isHeat(false);

  if (j.contains("boundaryConditions")) {
    for (const auto& bc : j["boundaryConditions"]) {
      BoundaryCondition boundaryCondition;
      if (bc.contains("name")) {
        boundaryCondition.name = bc["name"];
        isName = true;
      }
      if (bc.contains("type")) {
        boundaryCondition.type = bc["type"];

        if (bc["type"] == "Dirichlet") {
          isDirichlet = true;
        } else if (bc["type"] == "Neumann") {
          isNuemann = true;
        } else if (bc["type"] == "Third") {
          isThird = true;
        } else {
          throw std::runtime_error(
              "Error: type is not found. Dirichlet, Neumann, Third are "
              "available.");
        }
        isType = true;
      } else {
        throw std::runtime_error("Error: boundary type is not found.");
      }
      if (bc.contains("kind")) {
        boundaryCondition.kind = bc["kind"];
        isKind = true;
        if (bc["kind"] == "Heat") {
          isHeat = true;
        } else {
          throw std::runtime_error(
              "Error: kind is not found. Heat are "
              "available.");
        }
      } else {
        throw std::runtime_error("Error: boundary kind is not found.");
      }
      if (bc.contains("value")) {
        boundaryCondition.value = bc["value"].get<double>();
        isValue = true;
      }
      if (bc.contains("flux")) {
        boundaryCondition.flux = bc["flux"].get<double>();
        isFlux = true;
      }
      if (!isName) {
        throw std::runtime_error("Error: name is not found.");
      }
      if (!isKind) {
        throw std::runtime_error("Error: material kind is not found.");
      }
      if (!isType) {
        throw std::runtime_error("Error: type is not found.");
      }
      if (isDirichlet && !isValue) {
        throw std::runtime_error("Error: value is not found.");
      }
      if (isNuemann && !isFlux) {
        throw std::runtime_error("Error: flux is not found.");
      }
      if (isThird && !isValue) {
        throw std::runtime_error("Error: value is not found.");
      }
      if (isThird && !isFlux) {
        throw std::runtime_error("Error: flux is not found.");
      }
      if (isDirichlet && isNuemann) {
        throw std::runtime_error(
            "Error: Dirichlet and Neumann are defined at the same time.");
      }
      if (isDirichlet && isThird) {
        throw std::runtime_error(
            "Error: Dirichlet and Third are defined at the same time.");
      }
      if (isNuemann && isThird) {
        throw std::runtime_error(
            "Error: Neumann and Third are defined at the same time.");
      }
      if (isDirichlet && isNuemann && isThird) {
        throw std::runtime_error(
            "Error: Dirichlet, Neumann and Third are defined at the same "
            "time.");
      }
      boundaryConditions.push_back(boundaryCondition);
      isDirichlet = false;
      isNuemann = false;
      isThird = false;
      isName = false;
      isKind = false;
      isType = false;
      isFlux = false;
      isValue = false;
    }
  }

  if (j.contains("volumeProperties")) {
    for (const auto& vp : j["volumeProperties"]) {
      VolumeProperty volumeProperty;
      if (vp.contains("name")) {
        volumeProperty.name = vp["name"];
        isName = true;
      }
      if (vp.contains("kind")) {
        volumeProperty.kind = vp["kind"];
        isKind = true;
      }
      if (vp.contains("material")) {
        volumeProperty.material = vp["material"];
        isMaterial = true;
      }
      if (!isName) {
        throw std::runtime_error("Error: volume name is not found.");
      }
      if (!isKind) {
        throw std::runtime_error("Error: volume kind is not found.");
      }
      if (!isMaterial) {
        throw std::runtime_error("Error: volume material is not found.");
      }
      isName = false;
      isKind = false;
      isMaterial = false;
      volumeProperties.push_back(volumeProperty);
    }
  }
  if (j.contains("materialProperties")) {
    for (const auto& mp : j["materialProperties"]) {
      MaterialProperty material;
      if (mp.contains("name")) {
        material.name = mp["name"];
        isName = true;
      } else {
        throw std::runtime_error("Error: material name is not found.");
      }
      if (mp.contains("kind")) {
        material.type = mp["kind"];
        isKind = true;
      } else {
        throw std::runtime_error("Error: material kind is not found.");
      }
      if (mp.contains("density")) {
        material.density = mp["density"].get<double>();
      } else {
        throw std::runtime_error("Error: material density is not found.");
      }
      if (mp.contains("specificHeat")) {
        material.specificHeat = mp["specificHeat"].get<double>();
      } else {
        throw std::runtime_error("Error: material specificHeat is not found.");
      }
      if (mp.contains("thermalConductivity")) {
        material.thermalConductivity = mp["thermalConductivity"].get<double>();
      } else {
        throw std::runtime_error(
            "Error: material thermalConductivity is not found.");
      }
      if (material.type != "solid" && material.type != "fluid") {
        throw std::runtime_error("Error: type is not found.");
      }

      materials.push_back(material);
    }
  }
}
}  // namespace FemFy
#endif  // _CONFIG_H_
