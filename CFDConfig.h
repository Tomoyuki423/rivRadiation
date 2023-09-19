#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <vector>
#include <map>

using json = nlohmann::json;

struct GridSize {
    int x, y, z;
};

struct BoundaryCondition {
    std::string name;
    std::string type;
    std::string face;
    std::map<std::string, double> value;
};

class CFDConfig {
public:
    std::string solver;
    GridSize gridSize;
    double timeStep;
    int iterations;
    std::vector<BoundaryCondition> boundaryConditions;

    void loadFromFile(const std::string& filepath) {
        std::ifstream inputFile(filepath);
        if (!inputFile) {
            std::cerr << "Could not open file: " << filepath << std::endl;
            return;
        }

        json j;
        inputFile >> j;

        solver = j["solver"];
        gridSize.x = j["gridSize"]["x"];
        gridSize.y = j["gridSize"]["y"];
        gridSize.z = j["gridSize"]["z"];
        timeStep = j["timeStep"];
        iterations = j["iterations"];

        for (const auto& bc : j["boundaryConditions"]) {
            BoundaryCondition boundaryCondition;
            boundaryCondition.name = bc["name"];
            boundaryCondition.type = bc["type"];
            boundaryCondition.face = bc["face"];
            boundaryCondition.value = bc["value"].get<std::map<std::string, double>>();
            boundaryConditions.push_back(boundaryCondition);
        }
    }
};
