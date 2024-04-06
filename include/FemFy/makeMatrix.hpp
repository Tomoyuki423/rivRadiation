#ifndef ROM_MATRIX_h
#define ROM_MATRIX_h
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

typedef Eigen::Triplet<double> TP;
typedef Eigen::SparseMatrix<double> spMat;
/// param@dof_elm : the number of node in Element
namespace FemFy {
/// @brief Compute the generic matrix
template <typename T>
spMat computeGenericMatrix(
    Mesh& mesh, GaussPointRepository<double>& repository, int integrationOrder,
    std::function<double(const std::vector<double>&,
                         const std::vector<std::vector<double>>&,
                         const vec3d<double>&, double)>
        computeElementMatrixEntry,
    double propertyValue) {
  std::vector<double> shapeFunctions;
  std::vector<std::vector<double>> dNdrst;
  std::vector<std::vector<double>> jacobian;
  std::vector<vec3d<double>> node_coordinates;
  double weight;
  double detJ(0);
  int localMatrixSize;
  std::vector<TP> tripletVec;
  std::vector<std::vector<double>> localMatrix;

  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& it : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      GaussPointSet<double> gaussPointSet =
          repository.getGaussPointSet(it->getType(), integrationOrder);
      localMatrixSize = it->getNodeCount();
      localMatrix.resize(localMatrixSize,
                         std::vector<double>(localMatrixSize, 0.0));
      node_coordinates.resize(it->getNodeCount());

      for (int nodecnt = 0, size = it->getNodeCount(); nodecnt < size;
           nodecnt++) {
        node_coordinates[nodecnt] =
            mesh.nodes[it->getNodeIds(nodecnt)].getPos();
      }

      for (GaussPoint<double>& gp : gaussPointSet.points) {
        shapeFunctions = it->computeShapeFunctions(gp.coordinates.coord[0],
                                                   gp.coordinates.coord[1],
                                                   gp.coordinates.coord[2]);
        weight = gp.weight;
        dNdrst = it->computeShapeFunctionsDerivative(gp.coordinates.coord[0],
                                                     gp.coordinates.coord[1],
                                                     gp.coordinates.coord[2]);
        jacobian = calculateJacobianMatrix(node_coordinates, dNdrst);
        detJ = calculateDeterminant(jacobian);

        for (size_t i = 0; i < shapeFunctions.size(); i++) {
          for (size_t j = 0; j < shapeFunctions.size(); j++) {
            localMatrix[i][j] +=
                computeElementMatrixEntry(shapeFunctions, dNdrst,
                                          node_coordinates[j], detJ) *
                weight * propertyValue;
          }
        }
      }

      for (int i = 0; i < localMatrixSize; i++) {
        for (int j = 0; j < localMatrixSize; j++) {
          tripletVec.push_back(TP(it->getNodeIds(i) - 1, it->getNodeIds(j) - 1,
                                  localMatrix[i][j]));
        }
      }
    }
  }
  spMat M(mesh.getNumberOfNodes(), mesh.getNumberOfNodes());
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return M;
}

spMat computeMassMatrix(Mesh& mesh, CFDConfig& config,
                        GaussPointRepository<double>& repository,
                        int integrationOrder = 3) {
  std::vector<double> shapeFunctions;
  std::vector<std::vector<double>> dNdrst;
  std::vector<std::vector<double>> jacobian;
  std::vector<vec3d<double>> node_coordinates;
  double weight;
  double detJ(0);
  double rho = 1.;
  double cp = 1.;
  int localMatrixSize;
  std::vector<TP> tripletVec;
  std::vector<std::vector<double>> localMatrix;
  // will fix after implement GaussPointRepository
  // for (auto& configvp : config.volumeProperties) {
  //   for (auto& volume : mesh.physicalGrToElements.at(configvp.id))
  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& it : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      rho = config.materials[mesh.volumeForMaterial.at(volumeIndex).second]
                .density;
      cp = config.materials[mesh.volumeForMaterial.at(volumeIndex).second]
               .specificHeat;
      // for (size_t index = 0, size = mesh.elements.size(); index <
      // size;
      //      index++) {
      //   if (checkVolume(mesh, index)) {
      GaussPointSet<double> gaussPointSet =
          repository.getGaussPointSet(it->getType(), integrationOrder);
      localMatrixSize = it->getNodeCount();
      localMatrix.resize(localMatrixSize,
                         std::vector<double>(localMatrixSize, 0.0));
      node_coordinates.resize(it->getNodeCount());
      for (int nodecnt = 0, size = it->getNodeCount(); nodecnt < size;
           nodecnt++) {
        node_coordinates[nodecnt] =
            mesh.nodes[it->getNodeIds(nodecnt)].getPos();
      }

      for (GaussPoint<double>& gp : gaussPointSet.points) {
        shapeFunctions = it->computeShapeFunctions(gp.coordinates.coord[0],
                                                   gp.coordinates.coord[1],
                                                   gp.coordinates.coord[2]);
        weight = gp.weight;
        //   // write jacobian
        dNdrst = it->computeShapeFunctionsDerivative(gp.coordinates.coord[0],
                                                     gp.coordinates.coord[1],
                                                     gp.coordinates.coord[2]);
        jacobian = calculateJacobianMatrix(node_coordinates, dNdrst);
        detJ = calculateDeterminant(jacobian);
        for (size_t i = 0; i < shapeFunctions.size(); i++) {
          for (size_t j = 0; j < shapeFunctions.size(); j++) {
            localMatrix[i][j] += shapeFunctions[i] * shapeFunctions[j] *
                                 weight * rho * cp * detJ;
            // std::cout << "shapeFunction i(=" << i << ")" << shapeFunctions[i]
            //           << std::endl;
            // std::cout << "shapeFunction j(=" << j << ")" << shapeFunctions[j]
            //           << std::endl;
            // std::cout << "position" << gp.coordinates.coord[0] << ","
            //           << gp.coordinates.coord[1] << ","
            //           << gp.coordinates.coord[2] << std::endl;
            // std::cout << "weight" << weight << std::endl;
            // std::cout << "rho" << rho << std::endl;
            // std::cout << "cp" << cp << std::endl;
            // std::cout << "detJ" << detJ << std::endl;
          }
        }
      }
      for (int i = 0; i < localMatrixSize; i++) {
        /// sparse matrix is 0-indexed so we need to subtract 1
        for (int j = 0; j < localMatrixSize; j++) {
          // std::cout << "i:" << i << "j:" << j
          //           << "localMatrix[i][j]:" << localMatrix[i][j] <<
          //           std::endl;
          // if (it->getNodeIds(i) - 1 == 4 && it->getNodeIds(j) - 1 == 4) {
          //   std::cout << "i:" << i << "j:" << j
          //             << "localMatrix[i][j]:" << localMatrix[i][j] <<
          //             std::endl;
          // }
          tripletVec.push_back(TP(it->getNodeIds(i) - 1, it->getNodeIds(j) - 1,
                                  localMatrix[i][j]));
        }
      }
    }
  }
  spMat M(mesh.getNumberOfNodes(), mesh.getNumberOfNodes());
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return M;
}

spMat computeDiffusionMatrix(Mesh& mesh, CFDConfig& config,
                             GaussPointRepository<double>& repository,
                             int integrationOrder = 3) {
  std::vector<double> shapeFunctions;
  std::vector<std::vector<double>> dNdrst;
  std::vector<std::vector<double>> jacobian;
  std::vector<std::vector<double>> invjacobian;
  std::vector<std::vector<double>> invJxDN;
  std::vector<vec3d<double>> node_coordinates;
  double weight;
  double detJ(0);
  double kapp = 1;
  int localMatrixSize;
  std::vector<TP> tripletVec;
  std::vector<std::vector<double>> localMatrix;
  // will fix after implement GaussPointRepository
  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& it : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      // for (size_t index = 0, size = mesh.elements.size(); index < size;
      // index++) {
      //   if (checkVolume(mesh, index)) {

      kapp = config.materials[mesh.volumeForMaterial.at(volumeIndex).second]
                 .thermalConductivity;
      GaussPointSet<double> gaussPointSet =
          repository.getGaussPointSet(it->getType(), integrationOrder);
      localMatrixSize = it->getNodeCount();
      localMatrix.resize(localMatrixSize,
                         std::vector<double>(localMatrixSize, 0.0));
      node_coordinates.resize(it->getNodeCount());
      for (int nodecnt = 0, size = it->getNodeCount(); nodecnt < size;
           nodecnt++) {
        node_coordinates[nodecnt] =
            mesh.nodes[it->getNodeIds(nodecnt)].getPos();
      }
      for (GaussPoint<double>& gp : gaussPointSet.points) {
        shapeFunctions = it->computeShapeFunctions(gp.coordinates.coord[0],
                                                   gp.coordinates.coord[1],
                                                   gp.coordinates.coord[2]);
        weight = gp.weight;
        //   // write jacobian
        dNdrst = it->computeShapeFunctionsDerivative(gp.coordinates.coord[0],
                                                     gp.coordinates.coord[1],
                                                     gp.coordinates.coord[2]);

        jacobian = calculateJacobianMatrix(node_coordinates, dNdrst);
        detJ = calculateDeterminant(jacobian);
        if (detJ == 0) {
          continue;
        }
        invjacobian = calculateInverseMatrix(jacobian);
        invJxDN.resize(shapeFunctions.size(), std::vector<double>(3, 0.0));
        for (size_t i = 0, size = shapeFunctions.size(); i < size; i++) {
          for (size_t j = 0; j < 3; j++) {
            invJxDN[i][j] = invjacobian[j][0] * dNdrst[i][0] +
                            invjacobian[j][1] * dNdrst[i][1] +
                            invjacobian[j][2] * dNdrst[i][2];
          }
        }
        for (size_t i = 0; i < shapeFunctions.size(); i++) {
          for (size_t j = 0; j < shapeFunctions.size(); j++) {
            localMatrix[i][j] +=
                (invJxDN[i][0] * invJxDN[j][0] + invJxDN[i][1] * invJxDN[j][1] +
                 invJxDN[i][2] * invJxDN[j][2]) *
                weight * kapp * detJ;
          }
        }
      }
      for (int i = 0; i < localMatrixSize; i++) {
        /// sparse matrix is 0-indexed so we need to subtract 1
        for (int j = 0; j < localMatrixSize; j++) {
          // std::cout << "i:" << i << "j:" << j
          //           << "localMatrix[i][j]:" << localMatrix[i][j] <<
          //           std::endl;
          tripletVec.push_back(TP(it->getNodeIds(i) - 1, it->getNodeIds(j) - 1,
                                  localMatrix[i][j]));
        }
      }
    }
  }
  spMat M(mesh.getNumberOfNodes(), mesh.getNumberOfNodes());
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return M;
}

spMat computeConvectionMatrix(Mesh& mesh,
                              GaussPointRepository<double>& repository,
                              int integrationOrder = 3) {
  std::vector<double> shapeFunctions;
  std::vector<std::vector<double>> dNdrst;
  std::vector<std::vector<double>> jacobian;
  std::vector<vec3d<double>> node_coordinates;
  double weight;
  double detJ(0);
  // example
  vec3d<double> velocity(0., 0., 10.0);
  int localMatrixSize;
  std::vector<TP> tripletVec;
  std::vector<std::vector<double>> localMatrix;
  // will fix after implement GaussPointRepository
  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& it : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      // for (size_t index = 0, size = mesh.elements.size(); index < size;
      // index++) {
      //   if (checkVolume(mesh, index)) {
      GaussPointSet<double> gaussPointSet =
          repository.getGaussPointSet(it->getType(), integrationOrder);
      localMatrixSize = it->getNodeCount();
      localMatrix.resize(localMatrixSize,
                         std::vector<double>(localMatrixSize, 0.0));
      node_coordinates.resize(it->getNodeCount());

      for (GaussPoint<double>& gp : gaussPointSet.points) {
        shapeFunctions = it->computeShapeFunctions(gp.coordinates.coord[0],
                                                   gp.coordinates.coord[1],
                                                   gp.coordinates.coord[2]);
        weight = gp.weight;
        //   // write jacobian
        dNdrst = it->computeShapeFunctionsDerivative(gp.coordinates.coord[0],
                                                     gp.coordinates.coord[1],
                                                     gp.coordinates.coord[2]);
        for (int nodecnt = 0, size = it->getNodeCount(); nodecnt < size;
             nodecnt++) {
          node_coordinates[nodecnt] =
              mesh.nodes[it->getNodeIds(nodecnt)].getPos();
        }
        jacobian = calculateJacobianMatrix(node_coordinates, dNdrst);
        detJ = calculateDeterminant(jacobian);
        for (int i = 0; i < shapeFunctions.size(); i++) {
          for (int j = 0; j < shapeFunctions.size(); j++) {
            localMatrix[i][j] += (velocity.coord[0] * dNdrst[j][0] +
                                  velocity.coord[1] * dNdrst[j][1] +
                                  velocity.coord[2] * dNdrst[j][2]) *
                                 shapeFunctions[i] * weight * detJ;
          }
        }
      }
      for (int i = 0; i < localMatrixSize; i++) {
        /// sparse matrix is 0-indexed so we need to subtract 1
        for (int j = 0; j < localMatrixSize; j++) {
          tripletVec.push_back(TP(it->getNodeIds(i) - 1, it->getNodeIds(j) - 1,
                                  localMatrix[i][j]));
        }
      }
    }
  }
  spMat M(mesh.getNumberOfNodes(), mesh.getNumberOfNodes());
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return M;
}

// 対流方向の投影長を求める関数
template <typename T>
double computeConvectionLength(const std::vector<vec3d<T>>& node_coordinates,
                               const vec3d<T>& velocity) {
  double maxLength = 0.0;

  // 四面体要素の各辺ベクトルを計算し、速度ベクトルに投影
  for (size_t i = 0; i < node_coordinates.size(); ++i) {
    for (size_t j = i + 1; j < node_coordinates.size(); ++j) {
      vec3d edgeVector = node_coordinates[j] - node_coordinates[i];
      double projectionLength =
          fabs(edgeVector.dot(velocity)) / velocity.norm();
      if (projectionLength > maxLength) {
        maxLength = projectionLength;
      }
    }
  }
  return maxLength;
}

// SUPGパラメータの計算関数
template <typename T>
double computeSUPGParameter(double meshSize, const vec3d<T>& u, double nyu) {
  return std::min(meshSize / (2 * u.norm()),
                  1 / (u.norm() / meshSize + 2 * nyu / (meshSize * meshSize)));
}

spMat computeMassSUPGTerm(Mesh& mesh, GaussPointRepository<double>& repository,
                          int integrationOrder = 3) {
  // 方程式の要素行列へのSUPG追加項の組み込み関数
  // void addSUPGTerm(MatrixXd & elementMatrix, const VectorXd& u, double tau,
  //                  const MatrixXd& gradN, double weight) {
  std::vector<double> shapeFunctions;
  std::vector<std::vector<double>> dNdrst;
  std::vector<std::vector<double>> jacobian;
  std::vector<vec3d<double>> node_coordinates;
  // Example
  vec3d<double> velocity(0., 0., 10.0);
  double rho = 1.;
  double cp = 1.;
  double nyu = 1;
  double weight;
  double detJ(0);
  int localMatrixSize;
  std::vector<TP> tripletVec;
  std::vector<std::vector<double>> localMatrix;
  double meshSize;
  double tau;
  // will fix after implement GaussPointRepository
  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& it : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      // for (size_t index = 0, size = mesh.elements.size(); index < size;
      // index++) {
      //   if (checkVolume(mesh, index)) {
      GaussPointSet<double> gaussPointSet =
          repository.getGaussPointSet(it->getType(), integrationOrder);
      localMatrixSize = it->getNodeCount();
      localMatrix.resize(localMatrixSize,
                         std::vector<double>(localMatrixSize, 0.0));
      node_coordinates.resize(it->getNodeCount());
      for (int nodecnt = 0, size = it->getNodeCount(); nodecnt < size;
           nodecnt++) {
        node_coordinates[nodecnt] =
            mesh.nodes[it->getNodeIds(nodecnt)].getPos();
      }
      meshSize = computeConvectionLength(node_coordinates, velocity);
      tau = computeSUPGParameter(meshSize, velocity, nyu);
      for (GaussPoint<double>& gp : gaussPointSet.points) {
        shapeFunctions = it->computeShapeFunctions(gp.coordinates.coord[0],
                                                   gp.coordinates.coord[1],
                                                   gp.coordinates.coord[2]);
        weight = gp.weight;
        //   // write jacobian
        dNdrst = it->computeShapeFunctionsDerivative(gp.coordinates.coord[0],
                                                     gp.coordinates.coord[1],
                                                     gp.coordinates.coord[2]);

        jacobian = calculateJacobianMatrix(node_coordinates, dNdrst);
        detJ = calculateDeterminant(jacobian);
        size_t sizeNbNode = dNdrst.size();

        for (size_t i = 0; i < sizeNbNode; i++) {
          vec3d<double> gradNi = {dNdrst[i][0], dNdrst[i][1], dNdrst[i][2]};
          double dotProduct = gradNi.dot(velocity);
          for (size_t j = 0; j < sizeNbNode; j++) {
            localMatrix[i][j] +=
                rho * cp * tau * dotProduct * shapeFunctions[j] * weight * detJ;
          }
        }
      }
      for (int i = 0; i < localMatrixSize; i++) {
        /// sparse matrix is 0-indexed so we need to subtract 1
        for (int j = 0; j < localMatrixSize; j++) {
          tripletVec.push_back(TP(it->getNodeIds(i) - 1, it->getNodeIds(j) - 1,
                                  localMatrix[i][j]));
        }
      }
    }
  }
  spMat M(mesh.getNumberOfNodes(), mesh.getNumberOfNodes());
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return M;
}

spMat computeSUPGTerm(Mesh& mesh, GaussPointRepository<double>& repository,
                      int integrationOrder = 3) {
  // 方程式の要素行列へのSUPG追加項の組み込み関数
  // void addSUPGTerm(MatrixXd & elementMatrix, const VectorXd& u, double tau,
  //                  const MatrixXd& gradN, double weight) {
  std::vector<double> shapeFunctions;
  std::vector<std::vector<double>> dNdrst;
  std::vector<std::vector<double>> jacobian;
  std::vector<vec3d<double>> node_coordinates;
  // Example
  vec3d<double> velocity(0., 0., 10.0);
  double nyu = 1;
  double weight;
  double detJ(0);
  int localMatrixSize;
  std::vector<TP> tripletVec;
  std::vector<std::vector<double>> localMatrix;
  double meshSize;
  double tau;
  // will fix after implement GaussPointRepository
  for (size_t volumeIndex = 0; volumeIndex < mesh.volumeForMaterial.size();
       volumeIndex++) {
    for (auto& it : mesh.physicalGrToElements.at(
             mesh.volumeForMaterial.at(volumeIndex).first)) {
      // for (size_t index = 0, size = mesh.elements.size(); index < size;
      // index++) {
      //   if (checkVolume(mesh, index)) {
      GaussPointSet<double> gaussPointSet =
          repository.getGaussPointSet(it->getType(), integrationOrder);
      localMatrixSize = it->getNodeCount();
      localMatrix.resize(localMatrixSize,
                         std::vector<double>(localMatrixSize, 0.0));
      node_coordinates.resize(it->getNodeCount());
      for (int nodecnt = 0, size = it->getNodeCount(); nodecnt < size;
           nodecnt++) {
        node_coordinates[nodecnt] =
            mesh.nodes[it->getNodeIds(nodecnt)].getPos();
      }
      meshSize = computeConvectionLength(node_coordinates, velocity);
      tau = computeSUPGParameter(meshSize, velocity, nyu);
      for (GaussPoint<double>& gp : gaussPointSet.points) {
        shapeFunctions = it->computeShapeFunctions(gp.coordinates.coord[0],
                                                   gp.coordinates.coord[1],
                                                   gp.coordinates.coord[2]);
        weight = gp.weight;
        //   // write jacobian
        dNdrst = it->computeShapeFunctionsDerivative(gp.coordinates.coord[0],
                                                     gp.coordinates.coord[1],
                                                     gp.coordinates.coord[2]);

        jacobian = calculateJacobianMatrix(node_coordinates, dNdrst);
        detJ = calculateDeterminant(jacobian);
        size_t sizeNbNode = dNdrst.size();
        size_t sizeDimention = dNdrst[0].size();

        for (size_t i = 0; i < sizeNbNode; i++) {
          for (size_t j = 0; j < sizeNbNode; j++) {
            // SUPG追加項の計算
            vec3d<double> gradNi = {dNdrst[i][0], dNdrst[i][1], dNdrst[i][2]};
            vec3d<double> gradNj = {dNdrst[j][0], dNdrst[j][1], dNdrst[j][2]};
            double dotProduct = gradNi.dot(velocity) * gradNj.dot(velocity);
            localMatrix[i][j] += tau * dotProduct * weight * detJ;
            // localMatrix[i][j] += tau * velocity.coord[k] * dNdrst[i][k] *
            //                      velocity.coord[l] * dNdrst[j][l] *
            //                      weight * detJ;
          }
        }
      }
      for (int i = 0; i < localMatrixSize; i++) {
        /// sparse matrix is 0-indexed so we need to subtract 1
        for (int j = 0; j < localMatrixSize; j++) {
          tripletVec.push_back(TP(it->getNodeIds(i) - 1, it->getNodeIds(j) - 1,
                                  localMatrix[i][j]));
        }
      }
    }
  }
  spMat M(mesh.getNumberOfNodes(), mesh.getNumberOfNodes());
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  return M;
}

void applyDirichletBCforSPMAT(Eigen::SparseMatrix<double>& matrix,
                              const std::vector<std::pair<int, int>>& nodes) {
  typedef Eigen::Triplet<double> T;
  std::vector<T> triplets;
  std::set<int> rowsSet;  // 行インデックスの集合を直接作成

  // ノードのペアから行インデックスをセットに追加（0ベースインデックスのため-1）
  for (const auto& node : nodes) {
    int rowIndex = node.first - 1;  // 0ベースインデックスに変換
    rowsSet.insert(rowIndex);
  }

  // 既存の行列要素を保持（Dirichlet条件を適用する行以外）
  for (int k = 0; k < matrix.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it) {
      // Dirichlet条件を適用する行または列以外の要素を保持
      if (rowsSet.find(it.row()) == rowsSet.end()) {
        triplets.push_back(T(it.row(), it.col(), it.value()));
      }
    }
  }

  // Dirichlet条件の適用：指定された行の対角項を1に、それ以外を0に
  for (const int row : rowsSet) {
    triplets.push_back(T(row, row, 1.0));  // 対角項を1に
  }

  // 行列を再構築
  matrix.setZero();  // 既存の内容をクリア
  matrix.setFromTriplets(triplets.begin(), triplets.end());
}

/// @brief ベクトルBにDirichlet条件を適用する
/// @param B ベクトルB
/// @param rows Dirichlet条件を適用する行インデックス（1ベース）
/// @param values Dirichlet条件の値

void applyDirichletBCforVectorXd(Eigen::VectorXd& B,
                                 const std::vector<std::pair<int, int>>& rows,
                                 const CFDConfig& config) {
  std::unordered_map<int, double> boundaryIndex;
  for (auto& configbc : config.boundaryConditions) {
    if (configbc.kind == "Heat" && configbc.type == "Dirichlet") {
      boundaryIndex[configbc.id] = configbc.value;
    }
  }

  for (const auto& row : rows) {
    int rowIndex =
        row.first - 1;  // Dirichlet条件を適用する行インデックス（0ベース）
    int key = row.second;

    auto it = boundaryIndex.find(key);
    if (it == boundaryIndex.end()) {
      throw std::runtime_error("Key not found in boundaryIndex.");
    }

    double value = it->second;

    if (rowIndex >= 0 && rowIndex < B.size()) {
      B(rowIndex) = value;  // ベクトルBの対応する要素に値を設定
    } else {
      throw std::out_of_range("Row index is out of bounds of vector B.");
    }
  }
}

Eigen::VectorXd computeFvector(const Mesh& mesh, const CFDConfig& config,
                               GaussPointRepository<double>& repository,
                               int integrationOrder = 2) {
  unsigned int NbSurfaces = mesh.getNumberOfSurfaces(config);
  unsigned int NbNodes = mesh.getNumberOfNodes();
  Eigen::VectorXd F = Eigen::VectorXd::Zero(NbNodes);
  std::vector<double> shapeFunctions;
  std::vector<std::vector<double>> dNdrst;
  std::vector<std::vector<double>> jacobian;
  std::vector<vec3d<double>> node_coordinates;
  double weight;
  double detJ(0);
  double fluxContribution;
  for (auto& it : config.boundaryConditions) {
    if (it.kind == "Heat" && it.type == "Neumann") {
      for (auto& surface : mesh.physicalGrToElements.at(it.id)) {
        GaussPointSet<double> gaussPointSet =
            repository.getGaussPointSet(surface->getType(), integrationOrder);
        node_coordinates.resize(surface->getNodeCount());
        // for (int i = 0, size = surface->getNodeCount(); i < size; i++) {
        //   node_coordinates[i] =
        //   mesh.nodes.at(surface->getNodeIds(i)).getPos();
        // }
        // for (GaussPoint<double>& gp : gaussPointSet.points) {
        //   shapeFunctions = surface->computeShapeFunctions(
        //       gp.coordinates.coord[0], gp.coordinates.coord[1],
        //       gp.coordinates.coord[2]);
        //   weight = gp.weight;
        //   dNdrst = surface->computeShapeFunctionsDerivative(
        //       gp.coordinates.coord[0], gp.coordinates.coord[1],
        //       gp.coordinates.coord[2]);
        //   jacobian = calculateJacobianMatrix(node_coordinates, dNdrst);

        // detJ = calculateDeterminant(jacobian);
        weight = 1.0 / double(surface->getNodeCount());
        for (size_t i = 0, size = surface->getNodeCount(); i < size; i++) {
          // fluxContribution = it.flux * shapeFunctions[i] * weight * detJ;
          // fluxContribution = it.flux * shapeFunctions[i] * weight *
          //                    surface->getDimensionalProperties();
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
