#ifndef _OLDMATRIX_h
#define _OLDMATRIX_h

#include <FemFy/FemFy.hpp>
#include <FemFy/GaussPoints.hpp>
#include <FemFy/Simulation.hpp>
#include <FemFy/mesh.hpp>
#include <FemFy/readmsh.hpp>
#include <FemFy/utils.hpp>

///////////////////////////////////////////
//////////
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> spMat;
namespace FemFy {
/// param@dof_elm : the number of node in Element
// spMat makeM(int dof_elm, std::vector<std::vector<double>>& node,
//             std::vector<std::vector<unsigned int>>& volume,
//             std::vector<double>& rho, std::vector<double>& cp) {
spMat makeM(Mesh& mesh, CFDConfig& config) {
  // This program is based on tetra mesh.so dof_elm = 4
  int dof_elm = 4;
  unsigned int NbNodes = mesh.getNumberOfNodes();
  unsigned int NbVolumes = mesh.getNumberOfVolumes();
  std::vector<T> tripletVec;
  std::vector<std::vector<double>> Me(dof_elm, std::vector<double>(dof_elm));
  double totalvolume = 0.;
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  double elementvolume;
  double rho, cp;
  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& volume : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      rho = config.materials[mesh.volumeForMaterial.at(volumeIndex).second]
                .density;
      cp = config.materials[mesh.volumeForMaterial.at(volumeIndex).second]
               .specificHeat;
      // #x1 ~x3, y1 ~y3,z1~z3
      x1 = mesh.nodes[volume->getNodeIds(0)].getPos(0);
      x2 = mesh.nodes[volume->getNodeIds(1)].getPos(0);
      x3 = mesh.nodes[volume->getNodeIds(2)].getPos(0);
      x4 = mesh.nodes[volume->getNodeIds(3)].getPos(0);
      y1 = mesh.nodes[volume->getNodeIds(0)].getPos(1);
      y2 = mesh.nodes[volume->getNodeIds(1)].getPos(1);
      y3 = mesh.nodes[volume->getNodeIds(2)].getPos(1);
      y4 = mesh.nodes[volume->getNodeIds(3)].getPos(1);
      z1 = mesh.nodes[volume->getNodeIds(0)].getPos(2);
      z2 = mesh.nodes[volume->getNodeIds(1)].getPos(2);
      z3 = mesh.nodes[volume->getNodeIds(2)].getPos(2);
      z4 = mesh.nodes[volume->getNodeIds(3)].getPos(2);
      elementvolume = volume->getDimensionalProperties();
      // elementvolume =
      //     (-x1 * y2 * z3 + x1 * y2 * z4 + x1 * y3 * z2 - x1 * y3 * z4 -
      //      x1 * y4 * z2 + x1 * y4 * z3 + x2 * y1 * z3 - x2 * y1 * z4 -
      //      x2 * y3 * z1 + x2 * y3 * z4 + x2 * y4 * z1 - x2 * y4 * z3 -
      //      x3 * y1 * z2 + x3 * y1 * z4 + x3 * y2 * z1 - x3 * y2 * z4 -
      //      x3 * y4 * z1 + x3 * y4 * z2 + x4 * y1 * z2 - x4 * y1 * z3 -
      //      x4 * y2 * z1 + x4 * y2 * z3 + x4 * y3 * z1 - x4 * y3 * z2) /
      //     6.;

      if (elementvolume < 0) {
        std::cout << ("negative volume") << std::endl;
        std::cout << elementvolume << std::endl;
        elementvolume *= -1;
      } else if (elementvolume == 0) {
        std::cout << "0 volume" << std::endl;
      }

      for (int i = 0; i < dof_elm; i++) {
        for (int j = 0; j < dof_elm; j++) {
          if (i == j) {
            // mm[i][j] = rho[volume[e][j]] * cp[volume[e][j]] * elementvolume
            // / 10.;
            Me[i][j] = rho * cp * elementvolume / 10.;
          } else {
            // mm[i][j] = rho[volume[e][j]] * cp[volume[e][j]] * elementvolume
            // / 20.;
            Me[i][j] = rho * cp * elementvolume / 20.;
          }
          // Me[e][i][j] = mm[i][j];
        }
      }
      for (int i = 0; i < dof_elm; i++) {
        for (int j = 0; j < dof_elm; j++) {
          tripletVec.push_back(T(volume->getNodeIds(i) - 1,
                                 volume->getNodeIds(j) - 1, Me[i][j]));
        }
      }
    }
  }

  std::cout << "total volume, " << totalvolume << std::endl;

  spMat M(mesh.getNumberOfNodes(), mesh.getNumberOfNodes());
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  // std::cout << M << std::endl;
  return M;
}

spMat makeK(Mesh& mesh, CFDConfig& config) {
  int dof_elm = 4;
  unsigned int NbNodes = mesh.getNumberOfNodes();
  unsigned int NbVolumes = mesh.getNumberOfVolumes();
  std::vector<double> bb(dof_elm);
  std::vector<double> cc(dof_elm);
  std::vector<double> dd(dof_elm);
  std::vector<std::vector<double>> Be(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> Ce(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> De(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> kk(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> Ke(dof_elm, std::vector<double>(dof_elm));
  std::vector<T> tripletVec;
  spMat K(NbNodes, NbNodes);
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  double elementvolume;
  double kapper;

  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& volume : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      kapper = config.materials[mesh.volumeForMaterial.at(volumeIndex).second]
                   .thermalConductivity;
      // #x1 ~x3, y1 ~y3,z1~z3
      x1 = mesh.nodes[volume->getNodeIds(0)].getPos(0);
      x2 = mesh.nodes[volume->getNodeIds(1)].getPos(0);
      x3 = mesh.nodes[volume->getNodeIds(2)].getPos(0);
      x4 = mesh.nodes[volume->getNodeIds(3)].getPos(0);
      y1 = mesh.nodes[volume->getNodeIds(0)].getPos(1);
      y2 = mesh.nodes[volume->getNodeIds(1)].getPos(1);
      y3 = mesh.nodes[volume->getNodeIds(2)].getPos(1);
      y4 = mesh.nodes[volume->getNodeIds(3)].getPos(1);
      z1 = mesh.nodes[volume->getNodeIds(0)].getPos(2);
      z2 = mesh.nodes[volume->getNodeIds(1)].getPos(2);
      z3 = mesh.nodes[volume->getNodeIds(2)].getPos(2);
      z4 = mesh.nodes[volume->getNodeIds(3)].getPos(2);
      elementvolume = volume->getDimensionalProperties();
      bb[0] = -(y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3) /
              elementvolume / 6;
      bb[1] = -(-y1 * z3 + y1 * z4 + y3 * z1 - y3 * z4 - y4 * z1 + y4 * z3) /
              elementvolume / 6;
      bb[2] = -(y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2) /
              elementvolume / 6;
      bb[3] = -(-y1 * z2 + y1 * z3 + y2 * z1 - y2 * z3 - y3 * z1 + y3 * z2) /
              elementvolume / 6;
      cc[0] = -(-x2 * z3 + x2 * z4 + x3 * z2 - x3 * z4 - x4 * z2 + x4 * z3) /
              elementvolume / 6;
      cc[1] = -(x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3) /
              elementvolume / 6;
      cc[2] = -(-x1 * z2 + x1 * z4 + x2 * z1 - x2 * z4 - x4 * z1 + x4 * z2) /
              elementvolume / 6;
      cc[3] = -(x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) /
              elementvolume / 6;
      dd[0] = -(x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3) /
              elementvolume / 6;
      dd[1] = -(-x1 * y3 + x1 * y4 + x3 * y1 - x3 * y4 - x4 * y1 + x4 * y3) /
              elementvolume / 6;
      dd[2] = -(x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2) /
              elementvolume / 6;
      dd[3] = -(-x1 * y2 + x1 * y3 + x2 * y1 - x2 * y3 - x3 * y1 + x3 * y2) /
              elementvolume / 6;
      for (int i = 0; i < dof_elm; i++) {
        for (int j = 0; j < dof_elm; j++) {
          Be[i][j] = bb[i] * bb[j];
          Ce[i][j] = cc[i] * cc[j];
          De[i][j] = dd[i] * dd[j];
        }
      }
      for (int i = 0; i < dof_elm; i++) {
        for (int j = 0; j < dof_elm; j++) {
          Ke[i][j] = elementvolume * (kapper * Be[i][j] + kapper * Ce[i][j] +
                                      kapper * De[i][j]);
        }
      }
      for (int i = 0; i < dof_elm; i++) {
        for (int j = 0; j < dof_elm; j++) {
          tripletVec.push_back(T(volume->getNodeIds(i) - 1,
                                 volume->getNodeIds(j) - 1, Ke[i][j]));
        }
      }
    }
  }
  K.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return K;
}

spMat makeS(Mesh& mesh, CFDConfig& config) {
  vec3d<double> velocity(0., 0., 2.0);
  int dof_elm = 4;
  unsigned int NbNodes = mesh.getNumberOfNodes();
  unsigned int NbVolumes = mesh.getNumberOfVolumes();
  std::vector<double> bb(dof_elm);
  std::vector<double> cc(dof_elm);
  std::vector<double> dd(dof_elm);
  std::vector<std::vector<double>> Be(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> Ce(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> De(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> Se(dof_elm, std::vector<double>(dof_elm));
  std::vector<T> tripletVec;
  spMat S(NbNodes, NbNodes);
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  double elementvolume;
  double kapper;

  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& volume : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      // #x1 ~x3, y1 ~y3,z1~z3
      x1 = mesh.nodes[volume->getNodeIds(0)].getPos(0);
      x2 = mesh.nodes[volume->getNodeIds(1)].getPos(0);
      x3 = mesh.nodes[volume->getNodeIds(2)].getPos(0);
      x4 = mesh.nodes[volume->getNodeIds(3)].getPos(0);
      y1 = mesh.nodes[volume->getNodeIds(0)].getPos(1);
      y2 = mesh.nodes[volume->getNodeIds(1)].getPos(1);
      y3 = mesh.nodes[volume->getNodeIds(2)].getPos(1);
      y4 = mesh.nodes[volume->getNodeIds(3)].getPos(1);
      z1 = mesh.nodes[volume->getNodeIds(0)].getPos(2);
      z2 = mesh.nodes[volume->getNodeIds(1)].getPos(2);
      z3 = mesh.nodes[volume->getNodeIds(2)].getPos(2);
      z4 = mesh.nodes[volume->getNodeIds(3)].getPos(2);
      elementvolume = volume->getDimensionalProperties();
      // elementvolume =
      //     (-x1 * y2 * z3 + x1 * y2 * z4 + x1 * y3 * z2 - x1 * y3 * z4 -
      //      x1 * y4 * z2 + x1 * y4 * z3 + x2 * y1 * z3 - x2 * y1 * z4 -
      //      x2 * y3 * z1 + x2 * y3 * z4 + x2 * y4 * z1 - x2 * y4 * z3 -
      //      x3 * y1 * z2 + x3 * y1 * z4 + x3 * y2 * z1 - x3 * y2 * z4 -
      //      x3 * y4 * z1 + x3 * y4 * z2 + x4 * y1 * z2 - x4 * y1 * z3 -
      //      x4 * y2 * z1 + x4 * y2 * z3 + x4 * y3 * z1 - x4 * y3 * z2) /
      //     6.;
      bb[0] = -(y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3) /
              elementvolume / 6;
      bb[1] = -(-y1 * z3 + y1 * z4 + y3 * z1 - y3 * z4 - y4 * z1 + y4 * z3) /
              elementvolume / 6;
      bb[2] = -(y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2) /
              elementvolume / 6;
      bb[3] = -(-y1 * z2 + y1 * z3 + y2 * z1 - y2 * z3 - y3 * z1 + y3 * z2) /
              elementvolume / 6;
      cc[0] = -(-x2 * z3 + x2 * z4 + x3 * z2 - x3 * z4 - x4 * z2 + x4 * z3) /
              elementvolume / 6;
      cc[1] = -(x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3) /
              elementvolume / 6;
      cc[2] = -(-x1 * z2 + x1 * z4 + x2 * z1 - x2 * z4 - x4 * z1 + x4 * z2) /
              elementvolume / 6;
      cc[3] = -(x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) /
              elementvolume / 6;
      dd[0] = -(x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3) /
              elementvolume / 6;
      dd[1] = -(-x1 * y3 + x1 * y4 + x3 * y1 - x3 * y4 - x4 * y1 + x4 * y3) /
              elementvolume / 6;
      dd[2] = -(x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2) /
              elementvolume / 6;
      dd[3] = -(-x1 * y2 + x1 * y3 + x2 * y1 - x2 * y3 - x3 * y1 + x3 * y2) /
              elementvolume / 6;
      // for (int i = 0; i < dof_elm; i++) {
      //   for (int j = 0; j < dof_elm; j++) {
      //     Be[i][j] = bb[i] * bb[j];
      //     Ce[i][j] = cc[i] * cc[j];
      //     De[i][j] = dd[i] * dd[j];
      //   }
      // }
      for (int i = 0; i < dof_elm; i++) {
        for (int j = 0; j < dof_elm; j++) {
          Se[i][j] =
              elementvolume / 4. *
              (velocity[0] * bb[j] + velocity[1] * cc[j] + velocity[2] * dd[j]);
        }
      }
      for (int i = 0; i < dof_elm; i++) {
        for (int j = 0; j < dof_elm; j++) {
          tripletVec.push_back(T(volume->getNodeIds(i) - 1,
                                 volume->getNodeIds(j) - 1, Se[i][j]));
        }
      }
    }
  }
  S.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return S;
}

Eigen::VectorXd makeF(const Mesh& mesh, const CFDConfig& config) {
  unsigned int NbSurfaces = mesh.getNumberOfSurfaces(config);
  unsigned int NbNodes = mesh.getNumberOfNodes();
  Eigen::VectorXd F = Eigen::VectorXd::Zero(NbNodes);
  std::vector<double> shapeFunctions;
  std::vector<std::vector<double>> dNdrst;
  std::vector<std::vector<double>> jacobian;
  std::vector<vec3d<double>> node_coordinates;
  double weight;
  double fluxContribution;
  for (auto& it : config.boundaryConditions) {
    if (it.kind == "Heat" && it.type == "Neumann") {
      for (auto& surface : mesh.physicalGrToElements.at(it.id)) {
        weight = 1.0 / double(surface->getNodeCount());
        for (size_t i = 0, size = surface->getNodeCount(); i < size; i++) {
          fluxContribution =
              it.flux * weight * surface->getDimensionalProperties();
          F[surface->getNodeIds(i) - 1] += fluxContribution;
        }
        // }
      }
    }
  }
  return F;
}
}  // namespace FemFy

#endif
