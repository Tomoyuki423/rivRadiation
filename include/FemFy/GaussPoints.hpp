
/// @file GaussPoints class and functions for Gauss points.

#ifndef _GAUSSPOINTS_H_
#define _GAUSSPOINTS_H_

#include <FemFy/FemFy.hpp>
#include <FemFy/Vector3d.hpp>
#include <FemFy/mesh.hpp>

namespace FemFy {

/// Define a function to generate Gauss points for a given element type and
/// integration order.
template <typename T>
void generateGaussPointsForTETRAHEDRON(std::vector<T>& weight,
                                       std::vector<vec3d<T>>& pos,
                                       const int integrationOrder);
template <typename T>
void generateGaussPointsForTRIANGLE(std::vector<T>& weight,
                                    std::vector<vec3d<T>>& pos,
                                    const int integrationOrder);

/// @brief A class representing a Gauss point.
/// @tparam T The type of the coordinates of the Gauss point.

template <typename T>
class GaussPoint {
 public:
  T weight;
  vec3d<T> coordinates;
  GaussPoint() {}
  GaussPoint(T weight, vec3d<T> coordinates)
      : weight(weight), coordinates(coordinates) {}
};

template <typename T>
class GaussPointSet {
 public:
  std::vector<GaussPoint<T>> points;
  GaussPointSet() {}
  GaussPointSet(const std::vector<GaussPoint<T>>& points) : points(points) {}
  void addPoint(const GaussPoint<T>& point) { points.push_back(point); }
  const std::vector<GaussPoint<T>>& getPoints() const { return points; }
};

/// @brief A class representing a repository of Gauss point sets.
template <typename T>
class GaussPointRepository {
 public:
  // ガウス点セットを特定のキー（例:
  // 要素タイプと積分オーダーの組み合わせ）に関連付けて保存
  std::map<std::pair<GmshElementType, int>, GaussPointSet<T>> repository;

  GaussPointSet<T> gaussPointRepositorySet(const GmshElementType type,
                                           const int integrationOrder);

  void addGaussPointSet(const GmshElementType key, int integrationOrder) {
    if (!hasGaussPointSet(key, integrationOrder)) {
      std::pair<GmshElementType, int> repoKey =
          std::make_pair(key, integrationOrder);
      repository[repoKey] = gaussPointRepositorySet(key, integrationOrder);
    }
  }

  /// ガウス点セットがリポジトリに存在するかどうかを確認
  /// @param key ガウス点セットのキー
  /// @return ガウス点セットが存在する場合はtrue、それ以外はfalse
  bool hasGaussPointSet(const GmshElementType& key,
                        int integrationOrder) const {
    std::pair<GmshElementType, int> repoKey =
        std::make_pair(key, integrationOrder);
    return repository.count(repoKey) > 0;
  }
  // 特定のキーに対応するガウス点セットを取得
  const GaussPointSet<T>& getGaussPointSet(const GmshElementType& key,
                                           int integrationOrder) const {
    std::pair<GmshElementType, int> repoKey =
        std::make_pair(key, integrationOrder);
    return repository.at(repoKey);
  }
};

template <typename T>
GaussPointSet<T> GaussPointRepository<T>::gaussPointRepositorySet(
    const GmshElementType type, const int integrationOrder) {
  std::vector<T> weight;
  std::vector<vec3d<T>> pos;
  std::vector<GaussPoint<T>> gaussPoints;
  switch (static_cast<int>(type)) {
    case static_cast<int>(GmshElementType::POINT_1NODE):
      break;
    case static_cast<int>(GmshElementType::LINE_2NODE):
      weight = {1., 1.};
      pos = {{-1 / sqrt(3.), 0, 0}, {1 / sqrt(3.), 0, 0}};
      break;
    case static_cast<int>(GmshElementType::LINE_3NODE):
      break;
    case static_cast<int>(GmshElementType::TRIANGLE_3NODE):
      if (integrationOrder > 3 || integrationOrder < 1) {
        std::cerr << "Unsupported pair (TYPE, integration order) : "
                     "TRIANGLE_3NODE, "
                  << integrationOrder << ")" << std::endl;
      } else {
        generateGaussPointsForTRIANGLE(weight, pos, integrationOrder);
      }
      break;
    case static_cast<int>(GmshElementType::TRIANGLE_6NODE):
      if (integrationOrder > 3 || integrationOrder < 1) {
        std::cerr << "Unsupported pair (TYPE, integration order) : "
                     "TRIANGLE_6NODE, "
                  << integrationOrder << ")" << std::endl;
      } else {
        generateGaussPointsForTRIANGLE(weight, pos, integrationOrder);
      }
      break;
    case static_cast<int>(GmshElementType::TETRAHEDRON_4NODE):
      if (integrationOrder > 3 || integrationOrder < 1) {
        std::cerr << "Unsupported pair (TYPE, integration order) : "
                     "TETRAHEDRON_4NODE, "
                  << integrationOrder << ")" << std::endl;
      } else {
        generateGaussPointsForTETRAHEDRON(weight, pos, integrationOrder);
      }
      break;
    case static_cast<int>(GmshElementType::TETRAHEDRON_10NODE):
      if (integrationOrder > 3 || integrationOrder < 1) {
        throw std::runtime_error(
            "Unsupported pair (TYPE, integration order) : "
            "TETRAHEDRON_10NODE, " +
            integrationOrder);
      } else {
        generateGaussPointsForTETRAHEDRON(weight, pos, integrationOrder);
      }
      break;
    case static_cast<int>(GmshElementType::QUAD_4NODE):
      weight = {1, 1, 1, 1};
      pos = {{-1 / sqrt(3.), -1 / sqrt(3.), 0},
             {1 / sqrt(3.), -1 / sqrt(3.), 0},
             {1 / sqrt(3.), 1 / sqrt(3.), 0},
             {-1 / sqrt(3.), 1 / sqrt(3.), 0}};
      break;
    case static_cast<int>(GmshElementType::QUAD_8NODE):
      break;
    case static_cast<int>(GmshElementType::PRISM_6NODE):
      break;
    case static_cast<int>(GmshElementType::PRISM_15NODE):
      break;
    case static_cast<int>(GmshElementType::HEXAHEDRON_8NODE):
      break;
    case static_cast<int>(GmshElementType::HEXAHEDRON_20NODE):
      break;
    case static_cast<int>(GmshElementType::PYRAMID_5NODE):
      break;
    case static_cast<int>(GmshElementType::PYRAMID_13NODE):
      break;
    default:
      std::cerr << "GmshElementType not supported" << std::endl;
      exit(1);
  }
  for (size_t i = 0; i < weight.size(); i++) {
    gaussPoints.push_back({weight[i], pos[i]});
  }
  GaussPointSet<T> gaussPointSet(gaussPoints);
  return gaussPointSet;
}

template <typename T>
void generateGaussPointsForTRIANGLE(std::vector<T>& weight,
                                    std::vector<vec3d<T>>& pos,
                                    const int integrationOrder) {
  switch (integrationOrder) {
    case 1:
      weight = {1. / 3.};
      pos = {{1. / 3., 1. / 3., 0.}};
      break;
    case 2:
      weight = {1.0 / 6., 1.0 / 6., 1.0 / 6.};
      pos = {{1.0 / 6., 1.0 / 6, 0},
             {2.0 / 3., 1.0 / 6, 0},
             {1.0 / 6., 2.0 / 3., 0}};
      break;
    case 3:
      weight = {-27. / 96., 25. / 96., 25. / 96., 25. / 96.};
      pos = {{1.0 / 3.0, 1.0 / 3.0, 0.0},  // 中心点
             {0.2, 0.2, 0.0},
             {0.6, 0.2, 0.0},
             {0.2, 0.6, 0.0}};

      break;
    // 更に高次のルールを追加可能
    default:
      std::cerr << "Unsupported integration order:" << integrationOrder
                << " for LINE_2NODE. please 1 or 2 or 3." << std::endl;
      exit(1);
  }
}

template <typename T>
void generateGaussPointsForTETRAHEDRON(std::vector<T>& weight,
                                       std::vector<vec3d<T>>& pos,
                                       const int integrationOrder) {
  switch (integrationOrder) {
    case 1:
      weight = {1. / 6.};
      pos = {{1. / 4, 1. / 4., 1. / 4.}};
      break;
    case 2:
      weight = {1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0};
      pos = {{0.138196601125011, 0.138196601125011, 0.138196601125011},
             {0.585410196624969, 0.138196601125011, 0.138196601125011},
             {0.138196601125011, 0.585410196624969, 0.138196601125011},
             {0.138196601125011, 0.138196601125011, 0.585410196624969}};
      break;
    case 3:
      weight = {-0.8 / 6.0, 0.45 / 6.0, 0.45 / 6.0, 0.45 / 6.0, 0.45 / 6.0};
      pos = {{1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0},  // 中心点
             {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
             {1.0 / 2.0, 1.0 / 6.0, 1.0 / 6.0},
             {1.0 / 6.0, 1.0 / 2.0, 1.0 / 6.0},
             {1.0 / 6.0, 1.0 / 6.0, 1.0 / 2.0}};
      break;
    // 更に高次のルールを追加可能
    default:
      std::cerr << "Unsupported integration order:" << integrationOrder

                << " for TETRAHEDRON. please 1 or 2 or 3." << std::endl;
      exit(1);
  }
}
// operator<<をオーバーロードする非メンバ関数
template <typename T>
std::ostream& operator<<(std::ostream& os, const GaussPoint<T>& point) {
  os << "Weight: " << point.weight << ", Coordinates: ("
     << point.coordinates.coord[0] << ", " << point.coordinates.coord[1] << ", "
     << point.coordinates.coord[2] << ")";
  return os;
}
}  // namespace FemFy

#endif  // _GAUSSPOINTS_H_
