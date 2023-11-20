#ifndef ROM_READGMSHFILE_h
#define ROM_READGMSHFILE_h

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// void init(const std::string file_name);
// void getfileextention(const std::string file_name);
// void getdelimiter(const std::string file_name);
// std::vector<std::string> split(const std::string& input);
// void deletespace(std::string& buf);
// void skip(std::string& str);
// template <typename type>
// bool getonedata(std::vector<type>& information);
// template <typename type>
// void inputdata(const std::string t, std::vector<type>& information);

void deleteSpace(std::string& buf) {
  size_t pos;
  while ((pos = buf.find_first_of(" \t")) != std::string::npos) {
    buf.erase(pos, 1);
  }
}

std::vector<std::string> split(const std::string& input, const char dlim) {
  std::stringstream stream(input);
  std::string field;
  std::vector<std::string> result;
  while (getline(stream, field, dlim)) {
    result.push_back(field);
  }
  return result;
}

template <typename TYPE>
void inputData(const std::string T, const char dlim,
               std::vector<TYPE>& information) {
  std::vector<std::string> T1;
  std::stringstream ss;
  T1 = split(T, dlim);
  information.resize(T1.size());
  for (long long i = 0; i < T1.size(); i++) {
    ss << T1.at(i);
    ss >> information.at(i);
    ss.clear();
    ss.str("");
  }
}

void skip(std::string& str) {
  size_t found;
  size_t newfound = found = str.find_first_of("*");
  bool flag(false);
  while (found != std::string::npos) {
    newfound = str.find_first_of("*", found + 1);
    if (newfound == found + 1) {
      str.erase(found);
    }
    found = newfound;
  }
}

bool readGmshFile(std::string& filename, std::vector<std::vector<double>>& node,
                  std::vector<std::vector<unsigned int>>& volume,
                  std::vector<std::vector<unsigned int>>& surface,
                  std::vector<unsigned int>& volumeGrIndex,
                  std::vector<unsigned int>& surfaceGrIndex) {
  std::ifstream ifs;
  const char dlim = ',';
  unsigned int filedatasize(0);
  std::string str;
  std::vector<double> nodedata;
  std::vector<unsigned int> surfacedata;
  std::vector<unsigned int> volumedata;
  bool isNode(false);
  bool isLine(false);
  bool isSurface(false);
  bool isVolume(false);
  bool isBoundary(false);
  unsigned int preGrIndex;
  unsigned int numberofdata;
  std::string errmessage;
  ifs.open(filename, std::ios::in);
  if (ifs.fail()) {
    throw std::runtime_error("Fail to open " + filename);
  }

  while (getline(ifs, str)) {
    skip(str);
    if (str.empty()) continue;
    if (str.find("*ELEMENT") != std::string::npos) {
      if (isSurface) {
        numberofdata = filedatasize - preGrIndex;
        // std::cout << numberofdata << std::endl;
        surfaceGrIndex.push_back(surfaceGrIndex.back() + numberofdata);
      } else if (isVolume) {
        numberofdata = filedatasize - preGrIndex;
        // std::cout << numberofdata << std::endl;
        volumeGrIndex.push_back(volumeGrIndex.back() + numberofdata);
      }
    }
    filedatasize++;
    if (str.find("*NODE") != std::string::npos) {
      isNode = true;
      isLine = false;
      isSurface = false;
      isVolume = false;
      isBoundary = false;
      continue;
    }
    if (str.find("*ELEMENT") != std::string::npos) {
      if (str.find("Line") != std::string::npos) {
        isNode = false;
        isLine = true;
        isSurface = false;
        isVolume = false;
        isBoundary = false;
        continue;
      } else if (str.find("Surface") != std::string::npos) {
        isNode = false;
        isLine = false;
        isSurface = true;
        isVolume = false;
        isBoundary = false;
        preGrIndex = filedatasize;
        if (surfaceGrIndex.size() == 0) {
          surfaceGrIndex.push_back(0);
        }
        continue;
      } else if (str.find("Volume") != std::string::npos) {
        isNode = false;
        isLine = false;
        isSurface = false;
        isVolume = true;
        isBoundary = false;
        preGrIndex = filedatasize;
        if (volumeGrIndex.size() == 0) {
          volumeGrIndex.push_back(0);
        }
        continue;
      } else {
        errmessage = "\nCheck line " + std::to_string(filedatasize) + " in" +
                     filename + "\n-> " + str;
        throw std::runtime_error(errmessage);
      }
    }
    if (str.find("*BOUNDARY") != std::string::npos) {
      std::cout << "BOUNDARYH" << std::endl;
      isNode = false;
      isLine = false;
      isSurface = false;
      isVolume = false;
      isBoundary = true;
      continue;
    }

    if (isNode) {
      // std::cout << "node," << str << std::endl;
      inputData(str, dlim, nodedata);
      node.push_back({nodedata[1], nodedata[2], nodedata[3]});
    }
    if (isSurface) {
      inputData(str, dlim, surfacedata);
      surface.push_back(
          {surfacedata[1] - 1, surfacedata[2] - 1, surfacedata[3] - 1});
    }
    if (isVolume) {
      inputData(str, dlim, volumedata);
      volume.push_back({volumedata[1] - 1, volumedata[2] - 1, volumedata[3] - 1,
                        volumedata[4] - 1});
    }
    if (isBoundary) {
    }
    // deleteSpace(str);
    // if (information[0] == name) {
    //   for (int i = 0; i < 3; i++) {
    //     information.resize(0);
    //     for (long long j = 0, size = information.size(); j < size; j++) {
    //       data.at(5 * i + j) = std::stod(information.at(j));
    //     }
    //   }
    //   break;
    // }
  }
  if (surfaceGrIndex.back() == surface.size()) {
    surfaceGrIndex.pop_back();
  }
  if (volumeGrIndex.back() == volume.size()) {
    volumeGrIndex.pop_back();
  }
  ifs.close();
  return true;
}

#endif
