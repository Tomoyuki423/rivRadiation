#ifndef ROM_READGMSHFILE_h
#define ROM_READGMSHFILE_h

#include <stdexcept>
#include <string>
#include <vector>

#include "FemFy/mesh.hpp"
#include "ext/simple_vtk.hpp"

static const unsigned int DIMENTION = 3;
static const unsigned int MAXNUMBEROFVTKFILE = 100000;

enum class VTK_Cell_TYPE {
  DUMMY,
  VERTEX,
  POLY_VERREX,
  LINE,
  POLY_LINE,
  TRIANGLE,
  TRIANGLE_STRIP,
  POLYGON,
  PIXEL,
  QUAD,
  TETRA,
  VOXEL,
  HEXAHEDRON,
  WEDGE,
  PYRAMID
};

class OutputVTK {
 public:
  OutputVTK(const FemFy::Mesh& mesh, const double dt, const double T,
            const unsigned int outputIterationInterval,
            std::string filenameheader_)
      : NbNodes(mesh.getNumberOfNodes()),
        NbVolumes(mesh.getNumberOfVolumes()),
        filenameheader(filenameheader_) {
    exNbOutputVTK = static_cast<unsigned int>(
                        T / dt / static_cast<double>(outputIterationInterval)) +
                    1;
    setparam(mesh);
  }
  ~OutputVTK() {}
  void output(const std::vector<double>& phi, const FemFy::Mesh& mesh,
              const int outputcnt);
  void outputseries(const double dt,
                    const unsigned int outputIterationInterval);

 private:
  void setparam(const FemFy::Mesh& mesh);
  int numberofcomponent;
  const unsigned int NbNodes;
  const unsigned int NbVolumes;
  unsigned int exNbOutputVTK;
  std::vector<double> position;
  std::vector<unsigned int> offsets;
  std::vector<unsigned int> connectivity;
  std::vector<unsigned int> type;
  std::string outputdir;
  std::string filenameheader;
};

void OutputVTK::setparam(const FemFy::Mesh& mesh) {
  // if (cell_type == VTK_Cell_TYPE::VERTEX)
  //   numberofcomponent = 1;
  // else if (cell_type == VTK_Cell_TYPE::LINE)
  //   numberofcomponent = 2;
  // else if (cell_type == VTK_Cell_TYPE::TRIANGLE)
  //   numberofcomponent = 3;
  // else if (cell_type == VTK_Cell_TYPE::PIXEL ||
  //          cell_type == VTK_Cell_TYPE::QUAD ||
  //          cell_type == VTK_Cell_TYPE::QUAD ||
  //          cell_type == VTK_Cell_TYPE::TETRA)
  //   numberofcomponent = 4;
  // else if (cell_type == VTK_Cell_TYPE::VOXEL ||
  //          cell_type == VTK_Cell_TYPE::HEXAHEDRON)
  //   numberofcomponent = 8;
  // else if (cell_type == VTK_Cell_TYPE::WEDGE)
  //   numberofcomponent = 6;
  // else if (cell_type == VTK_Cell_TYPE::PYRAMID)
  //   numberofcomponent = 5;
  // if (mesh.Element->getType() == FemFy::GmshElementType::TETRAHEDRON_10NODE)
  // {
  //   numberofcomponent = 10;
  //   cell_type = VTK_Cell_TYPE::TETRA;
  // } else if (mesh.Element->getType() ==
  //            FemFy::GmshElementType::TETRAHEDRON_4NODE) {
  //   numberofcomponent = 4;
  //   cell_type = VTK_Cell_TYPE::TETRA;
  //   // } else if(mesh.getType() == FemFy::MeshType::TRIANGLE) {
  //   //   cell_type = VTK_Cell_TYPE::TRIANGLE;
  //   // } else if(mesh.getType() == FemFy::MeshType::QUAD) {
  //   //   cell_type = VTK_Cell_TYPE::QUAD;
  // } else {
  //   throw std::out_of_range("the mesh type is not supported");
  // }
  numberofcomponent = 4;
  VTK_Cell_TYPE cell_type;
  cell_type = VTK_Cell_TYPE::TETRA;

  if (exNbOutputVTK > MAXNUMBEROFVTKFILE) {
    throw std::out_of_range(
        "the number of file are needed less than MAX NUMBER OF VTK "
        "FILE(defined by header of output vtk");
  }

  position.resize(mesh.getNumberOfNodes() * 3);
  offsets.resize(NbVolumes);
  type.resize(NbVolumes);

  int aloc(0);
  for (size_t i = 1; i <= NbNodes; i++) {
    aloc = mesh.nodes.at(i).getId() - 1;
    position[aloc * 3] = mesh.nodes.at(i).getPos(0);
    position[aloc * 3 + 1] = mesh.nodes.at(i).getPos(1);
    if (DIMENTION == 2)
      position[aloc * 3 + 2] = 0;
    else if (DIMENTION == 3)
      position[aloc * 3 + 2] = mesh.nodes.at(i).getPos(2);
    else
      std::exit(1);
  }

  for (size_t i = 0; i < NbVolumes; i++) {
    offsets[i] = (i + 1) * numberofcomponent;
  }

  int cnt = 0;
  if (numberofcomponent != -1) {
    connectivity.resize(NbVolumes * numberofcomponent);
    // connectivity.resize(5);
    for (size_t i = 0; i < mesh.volumeForMaterial.size(); i++) {
      for (auto& it :
           mesh.physicalGrToElements.at(mesh.volumeForMaterial.at(i).first)) {
        for (size_t j = 0; j < it->getNodeCount(); j++) {
          connectivity.at(cnt * numberofcomponent + j) = it->getNodeIds(j) - 1;
        }
        cnt++;
      }
    }
    // for (int i = 0; i < NbVolumes; i++) {
    //   for (int j = 0; j < numberofcomponent; j++) {
    //     connectivity.at(i * numberofcomponent + j) = volume[i][j];
    //   }
    // }
  }

  for (size_t i = 0; i < NbVolumes; i++) {
    type[i] = static_cast<unsigned int>(cell_type);  // recomended cast
  }
}

void OutputVTK::output(const std::vector<double>& phi, const FemFy::Mesh& mesh,
                       const int outputcnt) {
  SimpleVTK gen;
  std::string outputvtk;
  std::string sdatanumber;
  unsigned int NbZero;
  unsigned int deg;
  std::vector<int> nodeId(NbNodes, 0.);
  std::vector<double> volumeId(NbVolumes, 0.);
  for (deg = 0; deg < log10(MAXNUMBEROFVTKFILE); deg++) {
    if (exNbOutputVTK - 1 < pow(10, deg)) break;
  }
  for (int i = 0; i < log10(MAXNUMBEROFVTKFILE); i++) {
    if (outputcnt < pow(10, i)) {
      NbZero = deg - i;
      if (outputcnt == 0) NbZero--;
      break;
    }
  }
  for (size_t i = 0; i < NbZero; i++) {
    sdatanumber += "0";
  }
  sdatanumber = sdatanumber + std::to_string(outputcnt);
  outputvtk = filenameheader + sdatanumber;

  //////////////////
  for (auto it : mesh.nodes) {
    nodeId[it.second.getId() - 1] = it.second.getId();
  }
  std::vector<double> ctemp(NbVolumes, 0.);
  // set cell data
  for (size_t i = 0; i < NbVolumes; i++) {
    for (int j = 0; j < numberofcomponent; j++) {
      ctemp[i] += phi[connectivity.at(i * numberofcomponent + j)];
    }
    ctemp[i] *= 0.10;
  }
  gen.init();
  gen.beginVTK("UnstructuredGrid");
  gen.beginContent();
  gen.beginPiece();

  gen.setNumberOfPoints(NbNodes);
  gen.setNumberOfCells(NbVolumes);

  gen.beginPoints();
  gen.beginDataArray();
  gen.setNumberType("Float32");
  gen.setNumberOfComponents(3);  // even if 2d, set 3 (all z axis value are 0)
  gen.setFormat("ascii");
  gen.addVector(position);
  gen.endDataArray();
  gen.endPoints();

  gen.beginCells();
  gen.beginDataArray("offsets", "Int32", "ascii");
  gen.addVector(offsets);
  gen.endDataArray();
  gen.beginDataArray("connectivity", "Int32", "ascii");
  gen.addVector(connectivity);
  gen.endDataArray();
  gen.beginDataArray("types", "Int32", "ascii");
  gen.addVector(type);
  gen.endDataArray();
  gen.endCells();

  gen.beginPointData();
  gen.setScalars("Temperature");
  gen.beginDataArray("Temperature", "Float32", "ascii");
  gen.addVector(phi);
  gen.endDataArray();
  // gen.endPointData();

  // gen.beginPointData();
  gen.setScalars("nodeId");
  gen.beginDataArray("nodeId", "Int32", "ascii");
  gen.addVector(nodeId);
  gen.endDataArray();
  gen.endPointData();

  gen.beginCellData();
  gen.setScalars("cellTemperature");
  gen.beginDataArray("cellTemperature", "Float32", "ascii");
  gen.addVector(ctemp);
  gen.endDataArray();
  gen.endCellData();

  gen.endPiece();
  gen.endContent();
  gen.endVTK();
  gen.generate(outputvtk);
}

void OutputVTK::outputseries(const double dt,
                             const unsigned int outputIterationInterval) {
  std::string outputFileName = outputdir + filenameheader + ".vtu.series";
  std::ofstream output(outputFileName);
  std::string filename;
  while (!output) {
    std::cout << "Unable to open " << outputFileName
              << ". so please input filename for massflow by yourself"
              << std::endl;
    std::cin >> outputFileName;
    output.close();
    outputFileName = outputdir + outputFileName + ".vtu.series";
    output.open(outputFileName);
  }
  std::string s = "{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n";
  std::string fileheader = "    { \"name\" : \"" + filenameheader;
  std::string timeheader = "";
  std::string filetime = ".vtu\", \"time\" : ";
  std::string filetail = " },\n";
  std::string endfile = "  ]\n}";
  std::string outputvtk;
  std::string sdatanumber;
  unsigned int NbZero;
  unsigned int deg;
  unsigned int itr;
  for (deg = 0; deg < log10(MAXNUMBEROFVTKFILE); deg++) {
    if (exNbOutputVTK - 1 < pow(10, deg)) break;
  }
  for (unsigned int itr = 0; itr < exNbOutputVTK; itr++) {
    sdatanumber = "";
    for (int j = 0; j < log10(MAXNUMBEROFVTKFILE); j++) {
      if (itr < pow(10, j)) {
        NbZero = deg - j;
        if (itr == 0) NbZero--;
        break;
      }
    }
    for (unsigned int i = 0; i < NbZero; i++) {
      sdatanumber += "0";
    }
    sdatanumber = sdatanumber + std::to_string(itr);
    outputvtk = filenameheader + sdatanumber;
    s += fileheader + sdatanumber + filetime +
         std::to_string(itr * dt * outputIterationInterval) + filetail;
  }
  s.pop_back();
  s.pop_back();
  s += "\n" + endfile;
  output << s;
  output.close();
}

#endif
