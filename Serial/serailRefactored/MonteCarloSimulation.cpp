// MonteCarloSimulation.tpp

#include "MonteCarloSimulation.h"
#include <cstdlib>
#include <fstream>


MonteCarloSimulation::MonteCarloSimulation(float interactionStrength, int latticeSize,  float Tolerance, float T_MIN, float T_MAX, float T_STEP):
     lattice(interactionStrength, latticeSize),  
        randVect(),
        energyResults(),
        magnetizationResults(),
        monteCarloStepsResults(),
        temperatures(),
        tolerance(Tolerance),  // Initialize tolerance with the specified value
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        L(latticeSize),
        N(latticeSize * latticeSize)
{   
    randVect = std::make_unique<std::vector<int> >(N);
    // Initialize randVect with createRandVect function
    createRandVector();
    // Initialize energyResults with std::make_unique       
    energyResults = std::make_unique<std::vector<float> >();
    // Initialize magnetizationResults with std::make_unique
    magnetizationResults = std::make_unique<std::vector<float> >();
    // Initialize monteCarloStepsResults with std::make_unique
    monteCarloStepsResults = std::make_unique<std::vector<int> >();
    // Initialize temperatures with std::make_unique
    temperatures = std::make_unique<std::vector<float> >();                     

}


void MonteCarloSimulation::simulatePhaseTransition() {
    int deltaE ;
    int deltaM ;
    int step;
    float m = 0;
    float mExact = 1;
    float T = T_MIN;

    float error = 0;
    std::array<float, 2> prob;
    
    while (T < T_MAX) {
        prob[0] = std::exp(-4 * lattice.getInteractionEnergy() / T);
        prob[1] = std::exp(-8 * lattice.getInteractionEnergy() / T);
        step = 0;
        mExact = mexact(T);
        m = static_cast<float>(lattice.getMagnetization()) / N;
        error = std::abs(std::abs(m) - mExact);
        while (error > tolerance ) {
            deltaE = 0;
            deltaM = 0; 
            createRandVector();
            simulateStep(prob, lattice.getLattice(), deltaM, deltaE);
            lattice.incrementMagnetization(deltaM);
            lattice.incrementEnergy(deltaE);
            m = static_cast<float>(lattice.getMagnetization()) / N;
            error = std::abs(std::abs(m) - mExact);
            step++;
        
        }

        temperatures->push_back(T);
        T += T_STEP;
        monteCarloStepsResults->push_back(step);
        energyResults->push_back(1);
        magnetizationResults->push_back(abs(m));
        step = 0;
        lattice.restoreRandomLattice();
    }
}



void MonteCarloSimulation::translateMatrix(std::vector<int>& inputMatrix) {
    std::unique_ptr<int[]> localCopy(new int[N]);
    std::memcpy(localCopy.get(), inputMatrix.data(), N * sizeof(int));

    for (int i = 0; i < N; ++i) {
        // Check if the new indices are within bounds
        if ((i + L) < N) {
            if ((i + 1) % L != 0)
                inputMatrix[i + L + 1] = localCopy[i];
            else
                inputMatrix[i + 1] = localCopy[i];
        } else if ((i + 1) % L != 0) {
            inputMatrix[i + L + 1 - N] = localCopy[i];
        } else {
            inputMatrix[0] = localCopy[i];
        }
    }
    
}

void MonteCarloSimulation::flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) {
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


void MonteCarloSimulation::simulateStep(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E) {
    for (unsigned long int i = 0; i < (N); i++) {
        int n = (*randVect)[i];
        if (n != -1) {
            flip(lattice, prob, n, M, E);
        }
    }
}


void MonteCarloSimulation::storeResultsToFile() const {
    // Open the file for writing
    std::ofstream outFile("result_" + std::to_string(N) + ".txt");

    // Check if the file is open
    if (!outFile.is_open()) {
        // Handle error: unable to open the file
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    // Write column headers
    outFile << "E  M  T N" << std::endl;

    // Determine the number of results to write
    std::size_t numResults = energyResults->size(); //they all have same lenght
    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << (*energyResults)[i] << " " << (*magnetizationResults)[i] << " " << (*temperatures)[i] << " " << (*monteCarloStepsResults)[i] << std::endl;

    }

    // Close the file
    outFile.close();
}


