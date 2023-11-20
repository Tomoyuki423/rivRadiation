#ifndef _MATRIX_h
#define _MATRIX_h

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> spMat;
/// param@dof_elm : the number of node in Element
spMat makeM(int dof_elm, std::vector<std::vector<double>>& node,
            std::vector<std::vector<unsigned int>>& volume,
            std::vector<double>& rho, std::vector<double>& cp) {
  // This program is based on tetra mesh.so dof_elm = 4
  unsigned int NbNodes = node.size();
  unsigned int NbVolumes = volume.size();
  std::vector<T> tripletVec;
  spMat M(NbNodes, NbNodes);
  std::vector<std::vector<double>> mm(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<std::vector<double>>> Me(
      NbVolumes,
      std::vector<std::vector<double>>(dof_elm, std::vector<double>(dof_elm)));
  double totalvolume = 0.;
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  double elementvolume;
  for (int e = 0; e < NbVolumes; e++) {
    // #x1 ~x3, y1 ~y3,z1~z3
    x1 = node[(volume[e][0])][0];
    x2 = node[(volume[e][1])][0];
    x3 = node[(volume[e][2])][0];
    x4 = node[(volume[e][3])][0];
    y1 = node[(volume[e][0])][1];
    y2 = node[(volume[e][1])][1];
    y3 = node[(volume[e][2])][1];
    y4 = node[(volume[e][3])][1];
    z1 = node[(volume[e][0])][2];
    z2 = node[(volume[e][1])][2];
    z3 = node[(volume[e][2])][2];
    z4 = node[(volume[e][3])][2];
    elementvolume =
        (-x1 * y2 * z3 + x1 * y2 * z4 + x1 * y3 * z2 - x1 * y3 * z4 -
         x1 * y4 * z2 + x1 * y4 * z3 + x2 * y1 * z3 - x2 * y1 * z4 -
         x2 * y3 * z1 + x2 * y3 * z4 + x2 * y4 * z1 - x2 * y4 * z3 -
         x3 * y1 * z2 + x3 * y1 * z4 + x3 * y2 * z1 - x3 * y2 * z4 -
         x3 * y4 * z1 + x3 * y4 * z2 + x4 * y1 * z2 - x4 * y1 * z3 -
         x4 * y2 * z1 + x4 * y2 * z3 + x4 * y3 * z1 - x4 * y3 * z2) /
        6.;

    if (elementvolume < 0) {
      std::cout << ("negative volume") << std::endl;
      std::cout << elementvolume << std::endl;
      elementvolume *= -1;
    } else if (elementvolume == 0) {
      std::cout << "0 volume" << std::endl;
    }

    totalvolume += elementvolume;
    for (int i = 0; i < dof_elm; i++) {
      for (int j = 0; j < dof_elm; j++) {
        if (i == j) {
          // mm[i][j] = rho[volume[e][j]] * cp[volume[e][j]] * elementvolume
          // / 10.;
          Me[e][i][j] = rho[e] * cp[e] * elementvolume / 10.;
        } else {
          // mm[i][j] = rho[volume[e][j]] * cp[volume[e][j]] * elementvolume
          // / 20.;
          Me[e][i][j] = rho[e] * cp[e] * elementvolume / 20.;
        }
        // Me[e][i][j] = mm[i][j];
      }
    }
  }

  std::cout << "total volume, " << totalvolume << std::endl;

  // #Setting matrix M
  for (int e = 0; e < NbVolumes; e++) {
    for (int i = 0; i < dof_elm; i++) {
      for (int j = 0; j < dof_elm; j++) {
        tripletVec.push_back(T(volume[e][i], volume[e][j], Me[e][i][j]));
      }
    }
  }
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  // std::cout << M << std::endl;
  return M;
}

spMat makeK(int dof_elm, std::vector<std::vector<double>>& node,
            std::vector<std::vector<unsigned int>>& volume,
            std::vector<double>& kapper) {
  unsigned int NbNodes = node.size();
  unsigned int NbVolumes = volume.size();
  std::vector<double> bb(dof_elm);
  std::vector<double> cc(dof_elm);
  std::vector<double> dd(dof_elm);
  std::vector<std::vector<double>> Be(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> Ce(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> De(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<double>> kk(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<std::vector<double>>> Ke(
      NbVolumes,
      std::vector<std::vector<double>>(dof_elm, std::vector<double>(dof_elm)));
  std::vector<T> tripletVec;
  spMat K(NbNodes, NbNodes);
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  double elementvolume;

  for (int e = 0; e < NbVolumes; e++) {
    // #x1 ~x3, y1 ~y3, z1~z3
    x1 = node[(volume[e][0])][0];
    x2 = node[(volume[e][1])][0];
    x3 = node[(volume[e][2])][0];
    x4 = node[(volume[e][3])][0];
    y1 = node[(volume[e][0])][1];
    y2 = node[(volume[e][1])][1];
    y3 = node[(volume[e][2])][1];
    y4 = node[(volume[e][3])][1];
    z1 = node[(volume[e][0])][2];
    z2 = node[(volume[e][1])][2];
    z3 = node[(volume[e][2])][2];
    z4 = node[(volume[e][3])][2];
    elementvolume =
        (-x1 * y2 * z3 + x1 * y2 * z4 + x1 * y3 * z2 - x1 * y3 * z4 -
         x1 * y4 * z2 + x1 * y4 * z3 + x2 * y1 * z3 - x2 * y1 * z4 -
         x2 * y3 * z1 + x2 * y3 * z4 + x2 * y4 * z1 - x2 * y4 * z3 -
         x3 * y1 * z2 + x3 * y1 * z4 + x3 * y2 * z1 - x3 * y2 * z4 -
         x3 * y4 * z1 + x3 * y4 * z2 + x4 * y1 * z2 - x4 * y1 * z3 -
         x4 * y2 * z1 + x4 * y2 * z3 + x4 * y3 * z1 - x4 * y3 * z2) /
        6.;
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
        Ke[e][i][j] =
            elementvolume * (kapper[e] * Be[i][j] + kapper[e] * Ce[i][j] +
                             kapper[e] * De[i][j]);
      }
    }
  }
  // #Ke
  // #Setting matrix K
  for (int e = 0; e < NbVolumes; e++) {
    for (int i = 0; i < dof_elm; i++) {
      for (int j = 0; j < dof_elm; j++) {
        tripletVec.push_back(T(volume[e][i], volume[e][j], Ke[e][i][j]));
      }
    }
  }
  K.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return K;
}

double calcTriangleArea(std::vector<double>& vec0, std::vector<double>& vec1,
                        std::vector<double>& vec2) {
  unsigned int dimention = vec1.size();
  double area;
  area = 0.5 *
         sqrt(pow(vec0[0] * vec1[1] - vec0[0] * vec2[1] - vec1[0] * vec0[1] +
                      vec1[0] * vec2[1] + vec2[0] * vec0[1] - vec2[0] * vec1[1],
                  2) +
              pow(vec0[0] * vec1[2] - vec0[0] * vec2[2] - vec1[0] * vec0[2] +
                      vec1[0] * vec2[2] + vec2[0] * vec0[2] - vec2[0] * vec1[2],
                  2) +
              pow(vec0[1] * vec1[2] - vec0[1] * vec2[2] - vec1[1] * vec0[2] +
                      vec1[1] * vec2[2] + vec2[1] * vec0[2] - vec2[1] * vec1[2],
                  2));
  return area;
}

std::vector<double> calSurfaceArea(
    int dof_elm, std::vector<std::vector<double>>& node,
    std::vector<std::vector<unsigned int>>& surface) {
  double totalArea = 0.;
  unsigned int NbNodes = node.size();
  unsigned int NbSurfaces = surface.size();
  std::vector<double> vec0(dof_elm, 0);
  std::vector<double> vec1(dof_elm, 0);
  std::vector<double> vec2(dof_elm, 0);
  std::vector<double> surfaceArea(NbSurfaces);

  for (int e = 0; e < NbSurfaces; e++) {
    for (int i = 0; i < dof_elm; i++) {
      vec0[i] = node[(surface[e][0])][i];
      vec1[i] = node[(surface[e][1])][i];
      vec2[i] = node[(surface[e][2])][i];
    }
    surfaceArea[e] = calcTriangleArea(vec0, vec1, vec2);
    totalArea += surfaceArea[e];
  }
  std::cout << "total area," << totalArea << std::endl;
  return surfaceArea;
}

std::vector<double> setSurfaceFlux(std::vector<unsigned int>& surfaceGrIndex,
                                   std::vector<double>& boundaryFluxCondition,
                                   std::vector<double> surfaceArea) {
  std::vector<double> surfaceFlux(surfaceArea.size());
  unsigned int NbSurfaces = surfaceArea.size();
  unsigned int NbSurfaceGrIndex = surfaceGrIndex.size();

  for (int e = 0; e < NbSurfaces; e++) {
    for (int i = NbSurfaceGrIndex - 1; i >= 0; i--) {
      if (e >= surfaceGrIndex[i]) {
        surfaceFlux[e] = boundaryFluxCondition[i];
        break;
      }
    }
  }
  return surfaceFlux;
}

Eigen::VectorXd makeF(int dof_elm, unsigned int NbNodes,
                      std::vector<std::vector<unsigned int>>& surface,
                      std::vector<double>& surfaceArea,
                      std::vector<double>& surfaceFlux) {
  unsigned int NbSurfaces = surface.size();
  Eigen::VectorXd F = Eigen::VectorXd::Zero(NbNodes);
  double weight = 1. / double(dof_elm);

  for (int e = 0; e < NbSurfaces; e++) {
    for (int i = 0; i < dof_elm; i++) {
      F[surface[e][i]] += weight * surfaceArea[e] * surfaceFlux[e];
    }
  }
  // std::cout << F << std::endl;
  return F;
}

#endif
