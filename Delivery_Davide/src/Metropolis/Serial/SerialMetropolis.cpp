
#include "SerialMetropolis.h"
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <random>


SerialMetropolis::SerialMetropolis(float interactionStrength, int latticeSize  , float T_MIN, float T_MAX, float T_STEP, long int IT)
: lattice(interactionStrength, latticeSize),  
        RandVect(),
        EnergyResults(),
        MagnetizationResults(),
        Temperatures(),
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        L(latticeSize),
        N(latticeSize * latticeSize),
        IT(IT)
      
{   

    // Initialize RandVect with std::make_unique
    RandVect = std::make_unique<std::vector<int> >();
    // Initialize rand_vect with create_rand_vect function
    create_rand_vector();
    // Initialize EnergyResults with std::make_unique       
    EnergyResults = std::make_unique<std::vector<float> >();
    // Initialize MagnetizationResults with std::make_unique
    MagnetizationResults = std::make_unique<std::vector<float> >();
    // Initialize Temperatures with std::make_unique
    Temperatures = std::make_unique<std::vector<float> >();


}


void SerialMetropolis::simulate_phase_transition() {
    int E_loc ;
    int M_loc ;
    int deltaE ;
    int deltaM ;
    int step;
    float m = 0;
    float T = T_MIN;
    int M = lattice.getMagnetization();
    float error = 0;
    std::array<float, 2> prob;

    while (T < T_MAX) {
        prob[0] = std::exp(-4 * lattice.getInteractionEnergy() / T);
        prob[1] = std::exp(-8 * lattice.getInteractionEnergy() / T);
        deltaE = 0;
        deltaM = 0;

        simulate_step(prob,lattice.get_lattice(), deltaM, deltaE);     
        lattice.incrementMagnetization(deltaM);
        lattice.incrementEnergy(deltaE);

        
        Temperatures->push_back(T);
        T += T_STEP;
        EnergyResults->push_back(1);
        m = static_cast<float>(lattice.getMagnetization()) / N;
        MagnetizationResults->push_back(abs(m));
        lattice.restore_random_lattice();
        }
    }   
    


void SerialMetropolis::create_rand_vector() {
    for (int j = 0;j<(IT);j++) {
        int r = rand() % L, c = rand() % L;
        RandVect->push_back ((r * L + c));
        }
}




void SerialMetropolis::flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) {
    int sum = 0;

    if (site < L) {
        sum += lattice[site + L * (L - 1)];
    } else {
        sum += lattice[site - L];
    }
    if (site % L == 0) {
        sum += lattice[site + (L - 1)];
    } else {
        sum += lattice[site - 1];
    }

    if (site >= L * (L - 1)) {
        sum += lattice[site - L * (L - 1)];
    } else {
        sum += lattice[site + L];
    }
    if ((site + 1) % L == 0) {
        sum += lattice[site - (L - 1)];
    } else {
        sum += lattice[site + 1];
    }
    int delta = 2 * sum * lattice[site];
    if (delta <= 0) {
        lattice[site] = -lattice[site];
    } else if (delta == 4) {
        float rnd = (rand() % 10000) / 1e4;
        if (rnd < prob[0]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    } else if (delta == 8) {
        float rnd = (rand() % 10000) / 1e4;
        if (rnd < prob[1]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    }
    M += 2 * lattice[site];
}


void SerialMetropolis::simulate_step (std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, int offset ) {
    for (unsigned long int i = 0; i < (IT);i++) {
        int n = (*RandVect)[i];
        flip(lattice, prob, n , M, E);
    }
}


void SerialMetropolis::store_results_to_file() const {

    std::string filePath = "results/result_" + std::to_string(N) + ".txt";
    // Open the file for writing
    std::ofstream outFile("result_" + std::to_string(N) + ".txt");

    // Check if the file is open
    if (!outFile.is_open()) {
        // Handle error: unable to open the file
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    // Write column headers
    outFile << "E  M  T " << std::endl;

    // Determine the number of results to write
    std::size_t numResults = EnergyResults->size(); //they all have same lenght
    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << (*EnergyResults)[i] << " " << (*MagnetizationResults)[i] << " " << (*Temperatures)[i]  << std::endl;

    }

    // Close the file
    outFile.close();
}


