#ifndef UTIL_H
#define UTIL_H

#include <FemFy/FemFy.hpp>
#include <FemFy/mesh.hpp>

/// utilties functions without class
namespace FemFy {

inline double calculateDeterminant2x2(
    const std::vector<std::vector<double>>& matrix) {
  return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}

inline double calculateDeterminant3x3(
    const std::vector<std::vector<double>>& matrix) {
  return matrix[0][0] *
             (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
         matrix[0][1] *
             (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
         matrix[0][2] *
             (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

inline double calculateDeterminant(
    const std::vector<std::vector<double>>& matrix) {
  int n = matrix.size();
  if (n == 1) {
    return matrix[0][0];
  } else if (n == 2) {
    return calculateDeterminant2x2(matrix);
  } else if (n == 3) {
    return calculateDeterminant3x3(matrix);
  }

  double determinant = 0;
  for (int i = 0; i < n; i++) {
    std::vector<std::vector<double>> minor(n - 1, std::vector<double>(n - 1));
    for (int j = 1; j < n; j++) {
      for (int k = 0; k < n; k++) {
        if (k < i) {
          minor[j - 1][k] = matrix[j][k];
        } else if (k > i) {
          minor[j - 1][k - 1] = matrix[j][k];
        }
      }
    }
    determinant +=
        (i % 2 == 0 ? 1 : -1) * matrix[0][i] * calculateDeterminant(minor);
  }

  return determinant;
}

inline std::vector<std::vector<double>> calculateJacobianMatrix(
    const std::vector<vec3d<double>>& node_coordinates,
    const std::vector<std::vector<double>>& dNdrst) {
  int n = dNdrst[0].size();  // dNdrstの列数に基づく
  std::vector<std::vector<double>> jacobian(
      n, std::vector<double>(n, 0.0));  // 3次元空間

  // Calculate the partial derivatives
  for (int i = 0; i < n; i++) {  // 形状関数の微分に対してループ
    for (size_t k = 0; k < node_coordinates.size();
         ++k) {                      // ノード座標に対してループ
      for (int j = 0; j < n; j++) {  // x, y, z 座標に対してループ
        jacobian[i][j] += dNdrst[k][i] * node_coordinates[k][j];
      }
    }
  }

  return jacobian;
}

// 3x3行列の逆行列を計算

inline std::vector<std::vector<double>> calculateInverseMatrix3x3(
    const std::vector<std::vector<double>>& matrix) {
  double det = calculateDeterminant(matrix);
  std::vector<std::vector<double>> inverse(3, std::vector<double>(3, 0.0));

  if (det == 0) {
    std::cout << "Matrix is singular, cannot calculate inverse." << std::endl;
    return inverse;  // 行列が特異であれば、逆行列は計算できません。
  }
  double invDet = 1.0 / det;
  inverse[0][0] =
      invDet * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]);
  inverse[0][1] =
      invDet * (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]);
  inverse[0][2] =
      invDet * (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]);
  inverse[1][0] =
      invDet * (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]);
  inverse[1][1] =
      invDet * (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]);
  inverse[1][2] =
      invDet * (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]);
  inverse[2][0] =
      invDet * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
  inverse[2][1] =
      invDet * (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]);
  inverse[2][2] =
      invDet * (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);

  return inverse;
}

inline std::vector<std::vector<double>> calculateInverseMatrix(
    const std::vector<std::vector<double>>& matrix) {
  if (matrix.size() == 3 && matrix[0].size() == 3) {
    return calculateInverseMatrix3x3(matrix);
  } else {
    throw std::invalid_argument(
        "Matrix size is not 3x3, cannot calculate inverse.");
  }
  return matrix;
}

}  // namespace FemFy
#endif  // UTIL_H
