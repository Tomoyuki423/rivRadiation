// Created by Tomoyuki Kawashima on 23/Sep/2023
#ifndef _FemFyROMDATA_H_
#define _FemFyROMATA_H_

#include <FemFy/FemFy.hpp>
#include <FemFy/makeMatrix.hpp>
#include <FemFy/mesh.hpp>
#include <FemFy/readmsh.hpp>
#include <FemFy/utils.hpp>

namespace FemFy {
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

std::vector<double> readAveData(const std::string filename) {
  std::string str;
  const char dlim = ',';
  std::vector<double> vecAvePhi;
  std::ifstream ifs;
  ifs.open(filename, std::ios::in);
  if (ifs.fail()) {
    throw std::runtime_error("Fail to open " + filename);
  }
  getline(ifs, str);
  inputData(str, dlim, vecAvePhi);
  ifs.close();
  return vecAvePhi;
}

std::vector<std::vector<double>> readSingularData(
    const std::string singularValueFile, const std::string singularVectorFile,
    const unsigned int NbNodes, const double threshold) {
  std::string str;
  const char dlim = ',';
  std::vector<double> singularvalue;
  std::string errmes;
  std::ifstream ifs;
  ifs.open(singularValueFile, std::ios::in);
  if (ifs.fail()) {
    throw std::runtime_error("Fail to open " + singularValueFile);
  }
  getline(ifs, str);
  inputData(str, dlim, singularvalue);
  ifs.close();
  double sum = std::accumulate(singularvalue.begin(), singularvalue.end(), 0.);
  double kiyoRate;
  unsigned int r = -1;
  for (int cnt = 0, size = singularvalue.size(); cnt < size; cnt++) {
    kiyoRate += singularvalue[cnt] / sum;
    if (kiyoRate > threshold) {
      r = cnt + 1;
      break;
    }
  }
  std::cout << "Select! the number of mode is " << r << std::endl;
  std::cout << "\nLoading singular value, now... " << std::endl;
  ifs.open(singularVectorFile, std::ios::in);
  if (ifs.fail()) {
    throw std::runtime_error("Fail to open " + singularVectorFile);
  }
  std::vector<std::vector<double>> vecMode(r, std::vector<double>(NbNodes));
  std::vector<double> modeData;
  unsigned int linecnt(0);
  while (getline(ifs, str)) {
    inputData(str, dlim, modeData);
    for (int i = 0; i < r; i++) {
      vecMode[i][linecnt] = modeData[i];
    }
    linecnt++;
  }
  ifs.close();
  return vecMode;
}

}  // namespace FemFy
#endif  // _FemFyROMDATA_H_
