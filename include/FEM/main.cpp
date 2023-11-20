//////////////////////////////////////////////////////////////////////
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

// #include <Eigen/src/SparseCholesky/SimplicialCholesky.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "./myFEM/makeMatrix.hpp"
#include "./myFEM/outputFile.hpp"
#include "./myFEM/readgmshfile.hpp"

std::chrono::system_clock::time_point ts_start, ts_end;
void start_eval(void);
void end_eval(void);
bool myFEM(std::string meshfile, std::string outputDir) {
  std::vector<std::vector<double>> node;
  std::vector<std::vector<unsigned int>> volume;
  std::vector<std::vector<unsigned int>> surface;
  std::vector<unsigned int> surfaceGrIndex;
  std::vector<unsigned int> volumeGrIndex;
  std::string errmesage;
  unsigned int NbNodes;
  unsigned int NbVolumes;
  unsigned int NbSurfaces;

  if (!readGmshFile(meshfile, node, volume, surface, volumeGrIndex,
                    surfaceGrIndex)) {
    errmesage = "Unable to open mesh file" + meshfile;
    throw std::runtime_error(errmesage);
  };

  // std::cout << "surface" << std::endl;
  // for (auto it : surfaceGrIndex) {
  //   std::cout << it << std::endl;
  // }
  // std::cout << "volume" << std::endl;
  // for (auto it : volumeGrIndex) {
  //   std::cout << it << std::endl;
  // }

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
  double dt = 0.01;
  double T = 20.;
  double t = 0;
  double theta = 0;
  unsigned int itr = 0;
  unsigned int outputcnt(0);
  unsigned int outputIterationInterval = 10;
  unsigned int outputConsoleInterval = 10;
  double outputTimeInterval = 10;
  std::string outputFileName;
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

  while (t < T) {
    preB = (M / dt - (1 - theta) * K);
    B = preB * phi;
    B = F + B;
    phi = chol.solve(B);
    itr++;
    t += dt;
    if (itr % outputConsoleInterval == 0) {
      if (itr % (10 * outputConsoleInterval) == 0) {
        std::cout << " Time, Temp.[K]" << std::endl;
      }
      std::cout << t << " ," << phi.mean() << std::endl;
    }
    if (itr % outputIterationInterval == 0) {
      outputcnt++;
      Eigen::Map<Eigen::VectorXd>(outputphi.data(), outputphi.size()) = phi;
      outputFileName =
          outputDir + "output" + std::to_string(outputcnt) + ".csv";
      outputResult(outputFileName, outputphi);
    }
  }
  return true;
}

int main() {
  try {
    std::string meshFile = "./msh/cube-Body_10_30.inp";
    std::string outputDir = "./result/";
    start_eval();
    bool flag = myFEM(meshFile, outputDir);
    end_eval();
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  return 0;
}

///
// Time calculation
///

void start_eval(void) { ts_start = std::chrono::system_clock::now(); }
void end_eval(void) {
  long sec, msec;
  ts_end = std::chrono::system_clock::now();  // ?ve??I??????
  sec = std::chrono::duration_cast<std::chrono::seconds>(ts_end - ts_start)
            .count();
  msec =
      std::chrono::duration_cast<std::chrono::milliseconds>(ts_end - ts_start)
          .count();
  if (sec < 0.01) {
    std::cout << msec << "[ms]" << std::endl;
  } else if (sec >= 0.01) {
    std::cout << sec << "[s]";
  }
}
