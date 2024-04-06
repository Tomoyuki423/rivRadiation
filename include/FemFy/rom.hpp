// Created by Tomoyuki Kawashima on 23/Sep/2023
#ifndef _FemFyROM_H_
#define _FemFyROM_H_

#include <FemFy/FemFy.hpp>
#include <FemFy/makeMatrix.hpp>
#include <FemFy/mesh.hpp>
#include <FemFy/readmsh.hpp>
#include <FemFy/utils.hpp>

namespace FemFy {
static const double THRESHOLD_SINGULAR_ACCUM = 0.999;
void end_eval(void);

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vector2matrix(
    const std::vector<std::vector<T>>& in) {
  unsigned int x = in.size();
  unsigned int y = in[0].size();
  std::vector<T> t;
  t.reserve(x * y);
  for (auto i : in) {
    t.insert(t.end(), i.begin(), i.end());
  }
  return (Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(&t[0], y,
                                                                       x))
      .transpose();
}

template <typename T>
Eigen::VectorXd romInitilize(
    Eigen::VectorXd& phi, Eigen::VectorXd& avePhi,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mode) {
  Eigen::VectorXd a(mode.size());
  a = mode * (phi - avePhi);
  return a;
}

bool myROM(std::string meshfile, std::string outputDir,
           std::string singularValueFile, std::string singularVectorFile,
           std::string aveValueFile, std::string outputVTKdir,
           std::string outputFilenameHeader) {
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

  M = makeM(4, node, volume, rho, cp);
  K = makeK(4, node, volume, Kapper);
  F = makeF(3, NbNodes, surface, surfaceArea, surfaceFlux);

  // get mode for ROM
  unsigned int NbMode(0);
  double threshold(THRESHOLD_SINGULAR_ACCUM);
  std::vector<std::vector<double>> vecMode;
  std::vector<double> vecAvePhi(NbNodes);
  vecAvePhi = readAveData(aveValueFile);
  vecMode = readSingularData(singularValueFile, singularVectorFile, NbNodes,
                             threshold);

  Eigen::VectorXd avePhi =
      Eigen::Map<Eigen::VectorXd>(vecAvePhi.data(), vecAvePhi.size());
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mode =
      vector2matrix(vecMode);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> romA;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> romVarK;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> romAveK;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> romB;
  Eigen::VectorXd romF;

  double dt = 0.01;
  double T = 1.0;
  double t = 0;

  std::string outputFileName;
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> preB;
  Eigen::VectorXd B(NbNodes);
  Eigen::VectorXd a(mode.size());
  A = M / dt + K * theta;
  romA = mode * A * mode.transpose();

  romAveK = mode * ((-(1 - theta) * K) * avePhi);
  romVarK = mode * (M / dt - (1 - theta) * K) * mode.transpose();
  romF = mode * F;
  // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
  Eigen::LDLT<Eigen::MatrixXd> chol;
  chol.compute(romA);
  phi.setConstant(100);
  a = romInitilize(phi, avePhi, mode);

  while (t < T) {
    romB = romF + romAveK + romVarK * a;
    a = chol.solve(romB);
    vtu.outputseries(dt, outputIterationInterval);
  }
  return true;
}

}  // namespace FemFy
#endif  // _FemFyROM_H_
