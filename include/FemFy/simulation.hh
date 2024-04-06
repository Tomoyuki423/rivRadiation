#ifndef _SIM_HH_
#define _SIM_HH_

#include <Eigen/src/Core/util/Constants.h>

#include <FemFy/Simulation.hpp>
#include <iostream>
#include <iterator>
#include <unordered_set>

namespace FemFy {

void Simulation::initialize() {
  std::cout << "set boundary condition id" << std::endl;
  setBoundaryConditionId();
  std::cout << "set volume property" << std::endl;
  setVolumeProperty();
  std::cout << "initialize gauss point" << std::endl;
  initializeGaussPointSet();
  std::cout << "initialize surface area" << std::endl;
  initializeSurfaceArea();
  std::cout << "initialize volume " << std::endl;
  initializeVolume();
}

void Simulation::initializeGaussPointSet() {
  // for (size_t i = 0, size = mesh.elements.size(); i < size; i++) {
  //   if (checkSurface(mesh, i)) {
  for (auto& configbc : config.boundaryConditions) {
    for (auto& surface : mesh.physicalGrToElements.at(configbc.id)) {
      integrationOrder = 2;
      // std::cout << "@simulation.hh in initializeGaussPointSet()" <<
      // std::endl; std::cout << "only surface integration order = " <<
      // integrationOrder
      //           << std::endl;
      repository.addGaussPointSet(surface->getType(), integrationOrder);
    }
  }
  for (auto& configvp : config.volumeProperties) {
    for (auto& volume : mesh.physicalGrToElements.at(configvp.id)) {
      // if (checkVolume(mesh, i)) {
      integrationOrder = 3;
      // std::cout << "@simulation.hh in initializeGaussPointSet()" <<
      // std::endl; std::cout << "only voulume integration order = " <<
      // integrationOrder
      //           << std::endl;
      repository.addGaussPointSet(volume->getType(), integrationOrder);
    }
  }
}

void Simulation::initializeSurfaceArea() {
  for (auto& configbc : config.boundaryConditions) {
    for (auto& surface : mesh.physicalGrToElements.at(configbc.id)) {
      surface->computeDimensionalProperties(mesh.nodes);
    }
    if (mesh.physicalGrMap.find("default") != mesh.physicalGrMap.end()) {
      for (auto& surface :
           mesh.physicalGrToElements.at(mesh.physicalGrMap.at("default"))) {
        surface->computeDimensionalProperties(mesh.nodes);
      }
    }
  }
  double sum = 0;
  int cnt = 0;
  for (auto& configbc : config.boundaryConditions) {
    for (auto& surface : mesh.physicalGrToElements.at(configbc.id)) {
      sum += surface->getDimensionalProperties();
      cnt++;
    }
    if (mesh.physicalGrMap.find("default") != mesh.physicalGrMap.end()) {
      for (auto& surface :
           mesh.physicalGrToElements.at(mesh.physicalGrMap.at("default"))) {
        sum += surface->getDimensionalProperties();
        cnt++;
      }
    }
  }
  std::cout << "@simulation.hh in initializeSurfaceArea()" << std::endl;
  std::cout << " Number of surface element, total surface area: " << cnt << ", "
            << sum << std::endl;
}

void Simulation::initializeVolume() {
  for (auto& configvp : config.volumeProperties) {
    for (auto& volume : mesh.physicalGrToElements.at(configvp.id)) {
      volume->computeDimensionalProperties(mesh.nodes);
    }
  }

  double sum = 0;
  int cnt = 0;
  for (auto& configvp : config.volumeProperties) {
    for (auto& volume : mesh.physicalGrToElements.at(configvp.id)) {
      sum += volume->getDimensionalProperties();
      cnt++;
    }
  }
  std::cout << "@simulation.hh in initializeVolume()" << std::endl;
  std::cout << "Number of volume element, total volume: " << cnt << "," << sum
            << std::endl;
}

template <typename DataType>
void Simulation::synchronizeGlobalVector(Eigen::VectorXd& globalValue,
                                         std::map<int, Node<DataType>>& nodes,
                                         DataType Node<DataType>::*attribute,
                                         bool toGlobal) {
  if (toGlobal) {
    for (const auto& node : nodes) {
      globalValue[node.second.getId() - 1] = node.second.*attribute;
    }
  } else {
    for (auto& node : nodes) {
      node.second.*attribute = globalValue[node.second.getId() - 1];
    }
  }
}

void Simulation::setVolumeProperty() {
  int i(0);
  for (auto& configmp : config.materials) {
    for (auto& configvp : config.volumeProperties) {
      if (configvp.material == configmp.name) {
        mesh.setVolumeForMaterial(configvp.name, i);
      }
      for (auto& it : mesh.physicalGrMap) {
        if (configvp.name == it.first) {
          configvp.id = it.second;
        }
      }
    }
    i++;
  }
}

void Simulation::setBoundaryConditionId() {
  for (auto& configbc : config.boundaryConditions) {
    for (auto& it : mesh.physicalGrMap) {
      if (configbc.name == it.first) {
        configbc.id = it.second;
      }
    }
  }
}

// std::vector<int> Simulation::setDirichletBcNode(const std::string kind) {
//   std::unordered_set<int> rowsSet;
//   std::vector<int> rows;
//   bool nameFound = false;

//   for (auto& configbc : config.boundaryConditions) {
//     if (configbc.kind == kind && configbc.type == "Dirichlet") {
//       for (auto it : mesh.physicalGrMap) {
//         if (configbc.name == it.first) {
//           nameFound = true;
//           int mapValue = it.second;  // mesh.physicalGrMap[it.first]の値
//           for (auto& surface :
//                mesh.physicalGrToElements[mesh.physicalGrMap[it.first]]) {
//             for (size_t i = 0, size = surface->getNodeCount(); i < size; i++)
//             {
//               rowsSet.insert(surface->getNodeIds(i));
//             }
//           }
//           break;
//         }
//       }
//       if (!nameFound) {
//         throw std::runtime_error(
//             "Error: configbc.name not found in mesh.physicalGrMap.");
//       }
//     }
//   }

//   // unordered_setからvectorに変換
//   rows.assign(rowsSet.begin(), rowsSet.end());
//   return rows;
// }
/// @brief setDirichletBcNode
/// @param kind
/// @return nodeId and physicalGrId.second

std::vector<std::pair<int, int>> Simulation::setDirichletBcNode(
    const std::string& kind) {
  std::unordered_set<int> rowsSet;
  std::vector<std::pair<int, int>> rows;
  bool nameFound = false;

  for (auto& configbc : config.boundaryConditions) {
    if (configbc.kind == kind && configbc.type == "Dirichlet") {
      for (auto& it : mesh.physicalGrMap) {
        if (configbc.name == it.first) {
          nameFound = true;
          int mapValue = it.second;  // mesh.physicalGrMap[it.first]の値
          for (auto& surface : mesh.physicalGrToElements[mapValue]) {
            for (size_t i = 0, size = surface->getNodeCount(); i < size; i++) {
              int nodeId = surface->getNodeIds(i);
              // もしrowsSetにnodeIdがまだ存在しないなら、rowsに追加
              if (rowsSet.insert(nodeId).second) {
                // std::cout << "nodeId: " << nodeId << ", mapValue: " <<
                // mapValue
                //           << std::endl;
                rows.push_back(std::make_pair(nodeId, mapValue));
              }
            }
          }
          break;
        }
      }
      if (!nameFound) {
        throw std::runtime_error(
            "Error: configbc.name not found in mesh.physicalGrMap.");
      }
    }
  }

  return rows;
}
}  // namespace FemFy

#endif  // _SIM_HH_
