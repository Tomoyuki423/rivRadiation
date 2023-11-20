// Created by Tomoyuki Kawashima on 23/Sep/2023
#ifndef _SIM_H_
#define _SIM_H_
#include "CFDConfig.h"
#include "element.h"

namespace FemFy {
class Simulation {
 private:
  CFDConfig config;
  void applyBoundaryConditions();  // 境界条件を適用
  void timeStep();                 // 1時間ステップを進める
  void outputResults();            // 結果を出力

  Mesh mesh_;              // メッシュオブジェクト
  Solver solver;           // ソルバーオブジェクト
  BoundaryConditions bc_;  // 境界条件オブジェクト
  Output output_;          // 出力オブジェクト
  double dt_;              // 時間ステップ
  int num_steps_;          // 総ステップ数
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

 public:
  Simulation() = default;
  // Simulation(const CFDConfig c, const std::vector<Node> n,
  //            const std::vector<Element*> elements)
  //     : config(c), nodes(n), elements(e) {}

  void initialize() {}
  void run() {
    // シミュレーションのロジックをここに書く
    // 設定は config メンバ変数からアクセスできる

    std::cout << "Running simulation with " << config.solver << " solver."
              << std::endl;
    std::cout << "Grid Size: (" << config.gridSize.x << ", "
              << config.gridSize.y << ", " << config.gridSize.z << ")"
              << std::endl;
    // ...
    solver.initialize();
  }
  void finalize() {}
};
}  // namespace FemFy
