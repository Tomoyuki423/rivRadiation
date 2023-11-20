#ifndef ROM_RESULT_H
#define ROM_RESULT_H

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void outputResult(const std::string outputFileName,
                  const std::vector<double> phi) {
  std::ofstream output(outputFileName);
  std::string filename;
  while (!output) {
    std::cout << "Unable to open " << outputFileName
              << ". so please input filename for massflow by yourself"
              << std::endl;
    std::cin >> filename;
    output.close();
    output.open(filename);
  }
  while (!Duct_Output) {
    std::cout << "Unable to open " << filename
              << ". so please input filename for massflow by yourself"
              << std::endl;
    std::cin >> filename;
    Duct_Output.close();
    Duct_Output.open(filename);
  }

  Duct_Output << "Duct name, flow direction, massflow[kg/s], Density[kg/m^3]"
              << std::endl;
  for (auto& it : duct) {
    Duct_Output << it.name_ << ',' << it.node_[0]->name_ << "->"
                << it.node_[1]->name_ << "," << it.massflow_ << ","
                << it.getRho() << std::endl;
  }
  Node_Output << "Node name" << std::flush;
  for (int i = 0; i < NUMBER_OF_SPECIES; i++) {
    Node_Output << "," << sp[i].getName() << std::flush;
  }
  Node_Output << " ,Enthalpy" << std::flush;
  Node_Output << " ,Temp.[K]" << std::flush;
  Node_Output << " ,Pressure.[Pa]" << std::flush;
  Node_Output << "\n" << std::flush;
  for (auto& it : node) {
    Node_Output << it.name_ << std::flush;
    for (int i = 0; i < NUMBER_OF_TRANSPORT; i++) {
      Node_Output << ", " << it.phi_[i] << std::flush;
    }
    Node_Output << ", " << it.temp_ << std::flush;
    Node_Output << ", " << it.pressure_ << std::flush;
    Node_Output << "\n" << std::flush;
  }
  Duct_Output.close();
  Node_Output.close();
  return true;
}
#endif  // ROM_RESULT_H
