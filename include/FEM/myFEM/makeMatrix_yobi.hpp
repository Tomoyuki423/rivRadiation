#ifndef ROM_MATRIX_h
#define ROM_MATRIX_h

#include <Eigen/Sparse>
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
          mm[i][j] = rho[volume[e][j]] * cp[volume[e][j]] * elementvolume / 10.;
        } else {
          mm[i][j] = rho[volume[e][j]] * cp[volume[e][j]] * elementvolume / 20.;
        }
        Me[e][i][j] = mm[i][j];
      }
    }
  }

  std::cout << "total volume " << totalvolume << std::endl;

  // #Setting matrix M
  for (int e = 0; e < NbVolumes; e++) {
    for (int i = 0; i < dof_elm; i++) {
      for (int j = 0; j < dof_elm; j++) {
        tripletVec.push_back(T(volume[e][i], volume[e][j], Me[e][i][j]));
      }
    }
  }
  M.setFromTriplets(tripletVec.begin(), tripletVec.end());
  std::cout << M << std::endl;
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
  std::vector<std::vector<double>> kk(dof_elm, std::vector<double>(dof_elm));
  std::vector<std::vector<std::vector<double>>> Ke(
      NbVolumes,
      std::vector<std::vector<double>>(dof_elm, std::vector<double>(dof_elm)));
  spMat K(dof_elm, dof_elm);
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  double elementvolume;

  for (int e = 0; i < NbVolumes; i++) {
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
    // Differential value of interpolating function
    // #x - direction is bb, y - direction is cc, z - direction is dd, \
    // N1 is[0], N2 is[1], N2 is[3], N3 is[4]
    bb[0] = (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    bb[1] = (-y1 * z3 + y1 * z4 + y3 * z1 - y3 * z4 - y4 * z1 + y4 * z3) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    bb[2] = (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    bb[3] = (-y1 * z2 + y1 * z3 + y2 * z1 - y2 * z3 - y3 * z1 + y3 * z2) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    cc[0] = (-x2 * z3 + x2 * z4 + x3 * z2 - x3 * z4 - x4 * z2 + x4 * z3) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    cc[1] = (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    cc[2] = (-x1 * z2 + x1 * z4 + x2 * z1 - x2 * z4 - x4 * z1 + x4 * z2) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    cc[3] = (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    dd[0] = (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    dd[1] = (-x1 * y3 + x1 * y4 + x3 * y1 - x3 * y4 - x4 * y1 + x4 * y3) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    dd[2] = (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    dd[3] = (-x1 * y2 + x1 * y3 + x2 * y1 - x2 * y3 - x3 * y1 + x3 * y2) /
            (x1 * y2 * z3 - x1 * y2 * z4 - x1 * y3 * z2 + x1 * y3 * z4 +
             x1 * y4 * z2 - x1 * y4 * z3 - x2 * y1 * z3 + x2 * y1 * z4 +
             x2 * y3 * z1 - x2 * y3 * z4 - x2 * y4 * z1 + x2 * y4 * z3 +
             x3 * y1 * z2 - x3 * y1 * z4 - x3 * y2 * z1 + x3 * y2 * z4 +
             x3 * y4 * z1 - x3 * y4 * z2 - x4 * y1 * z2 + x4 * y1 * z3 +
             x4 * y2 * z1 - x4 * y2 * z3 - x4 * y3 * z1 + x4 * y3 * z2);
    Be = np.outer(bb.T, bb);
    Ce = np.outer(cc.T, cc);
    De = np.outer(dd.T, dd);
    for (int i = 0; i < dof_elm; i++) {
      kk[:, i:i + 1] = volume * (k[element[e][i]] *
                                 Be[:, i:i + 1] + k[element[e][i]] * Ce
                                   [:, i:i + 1] + k[element[e][i]] * De
                                   [:, i:i + 1]);
    }
    if (Area == 0) {
      Area = 1e-6
    }
#print(Area)
    // Ke
    for i in range(dof_elm):
        for j in range(dof_elm):
#Ke[e][i][j] = kk[i][j]
        Ke[e]=kk
#Ke
    for e in range(NbElements):
        for i in range(dof_elm):
            for j in range(dof_elm):
                K[int(element[e][i])][int(element[e][j])]+=Ke[e][i][j]
    print("\n matrix K")
    print(K)
    return K
  }

//   makeF(dof_elm, node, boundaryFace, flux, boundaryFaceArea) {
//     double totalArea = 0.;
//     NbBoundaryFaces=len(boundaryFace)
//     NbNodes=len(node)
//     F=np.zeros(NbNodes)
//     weight  = 1./float(dof_elm)
// #FaceLabelArea = [0, 0, 0, 0, 0, 0]
//     x = np.empty(dof_elm)
//     y = np.empty(dof_elm)
//     z = np.empty(dof_elm)
//     if boundaryFaceArea[0]==-1:
//         for e in range(NbBoundaryFaces):
//             for i in range(dof_elm):
//                 x[i]=node[(boundaryFace[e][i])][0]
//                 y[i]=node[(boundaryFace[e][i])][1]
//                 z[i]=node[(boundaryFace[e][i])][2]
//             boundaryFaceArea[e] = calcArea(x,y,z)
//             totalArea+=boundaryFaceArea[e]
// #for i in range(6):
// #if boundaryFace[e][3] == i:
// #FaceLabelArea[i] += boundaryFaceArea[e]
// #break
//         print("total Area is " + str(totalArea))
// #print("Label Area is " + str(FaceLabelArea))

//     for e in range(NbBoundaryFaces):
//         for i in range(dof_elm):
//             F[boundaryFace[e][i]]+=weight*boundaryFaceArea[e]*flux[boundaryFace[e][3]]
// #bb, yy
// #print('Matrix F')
// #print(F)
//     return F,boundaryFaceArea

// double calcArea(vec1,vec2,vec3){
//       dimention =
//           len(x) if dimention == 3 : area =
//               0.5 * np.sqrt(np.abs(x[0] * y[1] - x[0] * y[2] - x[1] * y[0] +
//                                    x[1] * y[2] + x[2] * y[0] - x[2] * y[1]) *
//                                 *2 +
//                             np.abs(x[0] * z[1] - x[0] * z[2] - x[1] * z[0] +
//                                    x[1] * z[2] + x[2] * z[0] - x[2] * z[1]) *
//                                 *2 +
//                             np.abs(y[0] * z[1] - y[0] * z[2] - y[1] * z[0] +
//                                    y[1] * z[2] + y[2] * z[0] - y[2] * z[1]) *
//                                 *2)
// #print(str(e) + "," + str(boundaryFace[e][3]))
//                         return area
//     }
#endif
