#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream inputFile("tmp.dat");
    std::ofstream outputFile("filtered_E6.dat");
    
    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error opening files!" << std::endl;
        return 1;
    }
    
    std::string line;
    while (std::getline(inputFile, line)) {
        double x, y, z;
        if (std::sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z) == 3) {
            if (z == 1.89374466) {
                outputFile << line << std::endl;
            }
        }
    }
    
    inputFile.close();
    outputFile.close();
    
    std::cout << "Filtered data written to filtered_E6.dat" << std::endl;
    
    return 0;
}

