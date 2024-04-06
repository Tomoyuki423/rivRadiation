// Created by Tomoyuki Kawashima on 27/June/2023
// include guard
#ifndef SOLVER_H_
#define SOLVER_H_
#include <FemFy.hpp>
#include <vector3d.hpp>

namespace FemFy {

class Solver {
  Solver() = default;
  std::vector<std::vector<double>> node;
  std::vector<std::vector<unsigned int>> volume;
  std::vector<std::vector<unsigned int>> surface;
  std::vector<unsigned int> surfaceGrIndex;
  std::vector<unsigned int> volumeGrIndex;
  std::string errmesage;
  unsigned int NbNodes;
  unsigned int NbVolumes;
  unsigned int NbSurfaces;

  NbNodes = node.size();
  NbVolumes = volume.size();
  NbSurfaces = surface.size();
  std::cout << "Number of\nnodes, surfaces, volumes" << std::endl;
  std::cout << NbNodes << ", " << NbSurfaces << ", " << NbVolumes << std::endl;

  Eigen::SparseMatrix<double> M(NbNodes, NbNodes);
  Eigen::SparseMatrix<double> K(NbNodes, NbNodes);
  Eigen::VectorXd F(NbNodes);
  Eigen::VectorXd phi(NbNodes);
  std::vector<double> outputphi(NbNodes);

  // setting boundary condition
  std::vector<double> boundaryFluxCondition(surfaceGrIndex.size(), 0.);
  std::vector<double> surfaceFlux(NbSurfaces);
  std::vector<double> surfaceArea(NbSurfaces);

  //#Definition heat property
  std::vector<double> rho(NbVolumes, 1);
  std::vector<double> cp(NbVolumes, 1);
  std::vector<double> Kapper(NbVolumes, 1);

  // #Definition heat flux condition on boundary
  boundaryFluxCondition[4] = 10.;
  boundaryFluxCondition[5] = 5.;

  surfaceArea = calSurfaceArea(3, node, surface);
  surfaceFlux =
      setSurfaceFlux(surfaceGrIndex, boundaryFluxCondition, surfaceArea);

  M = makeM(4, node, volume, rho, cp);
  K = makeK(4, node, volume, Kapper);
  F = makeF(3, NbNodes, surface, surfaceArea, surfaceFlux);
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> preB;
  Eigen::VectorXd B(NbNodes);
  A = M / dt + K * theta;

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
  chol.compute(A);
  phi.setConstant(100);
  Eigen::Map<Eigen::VectorXd>(outputphi.data(), outputphi.size()) = phi;
  outputFileName = outputDir + "output" + std::to_string(outputcnt) + ".csv";
  outputResult(outputFileName, outputphi);
};
}  // namespace FemFy
