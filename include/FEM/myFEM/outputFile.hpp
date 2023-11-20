#ifndef ROM_RESULT_H
#define ROM_RESULT_H

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void outputResult(std::string outputFileName, const std::vector<double> phi) {
  std::ofstream output(outputFileName);
  std::string filename;
  while (!output) {
    std::cout << "Unable to open " << outputFileName
              << ". so please input filename for massflow by yourself"
              << std::endl;
    std::cin >> outputFileName;
    output.close();
    output.open(outputFileName);
  }

  for (auto& it : phi) {
    output << it << std::endl;
  }
  output.close();
}
#endif  // ROM_RESULT_H
