class Simulation {
private:
    CFDConfig config;

public:
    Simulation(const CFDConfig& c) : config(c) {}

    void run() {
        // シミュレーションのロジックをここに書く
        // 設定は config メンバ変数からアクセスできる

        std::cout << "Running simulation with " << config.solver << " solver." << std::endl;
        std::cout << "Grid Size: (" << config.gridSize.x << ", " << config.gridSize.y << ", " << config.gridSize.z << ")" << std::endl;
        // ...
    }
};
