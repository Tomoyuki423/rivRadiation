int main() {
    CFDConfig config;
    config.loadFromFile("config.json");

    Simulation sim(config);
    sim.run();

    return 0;
}
