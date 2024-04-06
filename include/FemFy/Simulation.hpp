// Created by Tomoyuki Kawashima on 23/Sep/2023
#ifndef _SIM_H_
#define _SIM_H_

#include <FemFy/FemFy.hpp>
#include <FemFy/makeMatrix.hpp>
#include <FemFy/mesh.hpp>
#include <FemFy/old_makeMatrix.hpp>
#include <FemFy/readmsh.hpp>
#include <FemFy/utils.hpp>

// #include "CFDConfig.h"

// #include <makeMatrix.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

// #include "GaussPoints.h"
#include "generate_vtu.hpp"
// #include "mesh.h"
// #include "utils.h"

typedef Eigen::SparseMatrix<double> spMat;
namespace FemFy {
void outputResult(std::string outputFileName, const std::vector<double> phi) {
  std::ofstream output(outputFileName);
  std::string filename;
  while (!output) {
    std::cout << "Unable to open " << outputFileName
              << ". so please input filename for massflow by yourself"
              << std::endl;
    std::cin >> outputFileName;
    output.close();
    output.open(outputFileName);
  }

  for (auto& it : phi) {
    output << it << std::endl;
  }
  output.close();
}

class Simulation {
 private:
  CFDConfig config;
  Mesh mesh;
  GaussPointRepository<double> repository;
  spMat massMatrix;
  spMat diffMatrix;
  spMat convectionMatrix;
  spMat supgMatrix;
  spMat supgMassMatrix;
  void timeStep();       // 1時間ステップを進める
  void outputResults();  // 結果を出力

  // Solver solver;          // ソルバーオブジェクト
  BoundaryCondition bc_;  // 境界条件オブジェクト
  // Output output_;         // 出力オブジェクト
  double dt_;      // 時間ステップ
  int num_steps_;  // 総ステップ数
  double dt = 0.01;
  double T = 20.;
  double t = 0;
  double theta = 1;
  unsigned int itr = 0;
  unsigned int outputcnt = 0;
  unsigned int outputIterationInterval = 10;
  unsigned int outputConsoleInterval = 10;
  double outputTimeInterval = 10;
  std::string outputFileName;
  std::vector<std::vector<std::pair<int, int>>> dirichletBCNodes;
  int integrationOrder;

 public:
  Simulation() = default;
  Simulation(const CFDConfig c, Mesh m) : config(c), mesh(std::move(m)) {
    config.gridSize.NbNode = mesh.nodes.size();
    config.gridSize.NbLine = 0;
    config.gridSize.NbSurface = 0;
    config.gridSize.NbVolume = 0;
  }

  void initialize();

  void initializeGaussPointSet();

  void initializeSurfaceArea();
  void initializeVolume();
  void setVolumeProperty();
  void setBoundaryConditionId();
  template <typename T>
  void synchronizeGlobalVector(Eigen::VectorXd& globalValue,
                               std::map<int, Node<T>>& nodes,
                               T Node<T>::*attribute, bool toGlobal);

  std::vector<std::pair<int, int>> setDirichletBcNode(const std::string& kind);

  void run() {
    int cntrun = 0;
    initialize();
    std::cout << "initialize done" << std::endl;
    integrationOrder = 3;
    std::cout << "check mass" << std::endl;
    // massMatrix = computeMassMatrix(mesh, config, repository,
    // integrationOrder);
    massMatrix = makeM(mesh, config);
    // diffMatrix =
    //     computeDiffusionMatrix(mesh, config, repository, integrationOrder);
    std::cout << "check diff" << std::endl;
    diffMatrix = makeK(mesh, config);
    // convectionMatrix =
    //     computeConvectionMatrix(mesh, repository, integrationOrder);
    convectionMatrix = makeS(mesh, config);
    // supgMatrix = computeSUPGTerm(mesh, repository, integrationOrder);
    // supgMassMatrix = computeMassSUPGTerm(mesh, repository, integrationOrder);
    Eigen::VectorXd F(mesh.getNumberOfNodes());
    F.setZero();
    Eigen::VectorXd phi(mesh.getNumberOfNodes());
    phi.setZero();
    std::vector<double> outputphi(mesh.getNumberOfNodes());
    Eigen::SparseMatrix<double> A;
    Eigen::SparseMatrix<double> preB;
    Eigen::VectorXd B(mesh.getNumberOfNodes());
    if (theta == 0) {
      A = massMatrix;
      // A = (massMatrix + supgMassMatrix);
    } else {
      // A = (massMatrix) + (diffMatrix)*dt * theta;
      A = (massMatrix) + (diffMatrix + convectionMatrix) * dt * theta;
      // A = (massMatrix + supgMassMatrix) + (diffMatrix + convectionMatrix +
      // supgMatrix) * dt * theta;
    }
    std::cout << "A:" << std::endl;
    std::cout << A << std::endl;
    std::cout << "M:" << std::endl;
    std::cout << massMatrix << std::endl;
    // std::cout << supgMassMatrix << std::endl;
    // std::cout << supgMatrix << std::endl;
    dirichletBCNodes.resize(1);
    dirichletBCNodes[0] = setDirichletBcNode("Heat");
    applyDirichletBCforSPMAT(A, dirichletBCNodes[0]);

    ///////////////////debug//////////////////////////
    // std::cout << A << std::endl;
    // std::cout << "@Simulation.hpp in run(), display node pos" << std::endl;
    // for (auto it : mesh.nodes) {
    //   std::cout << it.second.getPos(0) << "," << it.second.getPos(1) << ","
    //             << it.second.getPos(2) << std::endl;
    // }
    // for (auto& it : mesh.elements) {
    //   std::cout << "element id:" << it->getId() << std::endl;
    //   for (size_t i = 0; i < it->getNodeCount(); i++) {
    //     std::cout << "node id:" << it->getNodeIds(i) << std::endl;
    //   }
    // }

    // for (auto& it : dirichletBCNodes[0]) {
    //   std::cout << "check row :" << it << std::endl;
    // }
    //////////////////////////////////////////////////////

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> bicg;
    bicg.compute(A);
    phi.setConstant(0);
    for (auto& it : mesh.nodes) {
      // if (it.second.getPos(2) > 13.0 && it.second.getPos(2) < 17.0) {
      // if (it.second.getPos(2) < 15.0) {
      //   phi(it.second.getId() - 1) = 10;
      // }
    }

    Eigen::Map<Eigen::VectorXd>(outputphi.data(), outputphi.size()) = phi;
    outputFileName = "./output" + std::to_string(outputcnt) + ".csv";
    // outputResult(outputFileName, outputphi);
    std::cout << "check VTK:" << std::endl;
    OutputVTK outputVTK(mesh, dt, T, outputIterationInterval, "./output");
    outputVTK.output(outputphi, mesh, outputcnt);
    if (theta == 1) {
      preB = massMatrix;
    } else {
      // preB = (massMatrix) - (1 - theta) * dt * diffMatrix;
      preB = (massMatrix) - (1 - theta) * dt * (diffMatrix + convectionMatrix);
      // preB = (massMatrix + supgMassMatrix) -
      //        (1 - theta) * dt * (diffMatrix + convectionMatrix + supgMatrix);
    }

    // int boundaryintegrationOrder = 2;
    // F = computeFvector(mesh, config, repository, boundaryintegrationOrder);
    F = makeF(mesh, config);
    std::cout << "check before run:" << std::endl;
    while (t < T) {
      B = preB * phi;
      B = F * dt + B;
      applyDirichletBCforVectorXd(B, dirichletBCNodes[0], config);
      phi = bicg.solveWithGuess(B, phi);
      // synchronizeGlobalVector(phi, mesh.nodes, &Node<double>::temperature,
      //                         false);
      itr++;
      t += dt;
      if (itr % outputConsoleInterval == 0) {
        if (itr % (outputConsoleInterval) == 0) {
          std::cout << " Time, Temp.[K]" << std::endl;
        }
        std::cout << t << " ," << phi.mean() << std::endl;
      }
      if (itr % outputIterationInterval == 0) {
        outputcnt++;
        Eigen::Map<Eigen::VectorXd>(outputphi.data(), outputphi.size()) = phi;
        outputFileName = "./output" + std::to_string(outputcnt) + ".csv";
        // outputResult(outputFileName, outputphi);
        outputVTK.output(outputphi, mesh, outputcnt);
      }
    }
    // std::cout << massMatrix << std::endl;
    // std::cout << diffMatrix << std::endl;
    // std::cout << A << std::endl;
    // std::cout << B << std::endl;
    // std::cout << phi << std::endl;
    // std::cout << F << std::endl;
    // temporaly check
    std::cout << "physicalGrToElemetns:" << mesh.physicalGrToElements.size()
              << std::endl;
    std::cout << "nodes size:" << mesh.getNumberOfNodes() << std::endl;
    std::cout << "Elementss size:" << mesh.getNumberOfElements() << std::endl;
    std::cout << "Volume size:" << mesh.getNumberOfVolumes() << std::endl;
  }
  void finalize() {}
};
}  // namespace FemFy

#include <FemFy/simulation.hh>
#endif  // _SIM_H_
