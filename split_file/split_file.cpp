#include "../Solution_Items/GATE_data_structure.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> // For std::atoi


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <filename> <K>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    int K = std::atoi(argv[2]);
    
    if (K <= 0) {
        std::cerr << "Invalid value for K. It should be a positive integer." << std::endl;
        return 1;
    }

    std::ifstream inputFile(filename, std::ios::binary);
    std::ofstream outputFile("output.bin", std::ios::binary);

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the binary file: " << filename << std::endl;
        return 1;
    }

    if (!outputFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return 1;
    }

    mini_coinc_event data;
    int counter = 0;

    while (inputFile.read(reinterpret_cast<char*>(&data), sizeof(mini_coinc_event))) {
        counter++;
        if (counter % K == 0) {
            // Write the Kth record to the output file
            outputFile.write(reinterpret_cast<const char*>(&data), sizeof(mini_coinc_event));
        }
    }

    inputFile.close();
    outputFile.close();

    return 0;
}