// MonteCarloSimulation.tpp

#include "MonteCarloSimulation.h"
#include <cstdlib>
#include <fstream>



template <std::size_t N>
MonteCarloSimulation<N>::MonteCarloSimulation(float interactionStrength, int latticeSize,  float Tolerance, float T_MIN, float T_MAX, float T_STEP)
    : lattice(interactionStrength, latticeSize),  
        randVect(),
        energyResults(),
        magnetizationResults(),
        monteCarloStepsResults(),
        temperatures(),
        tolerance(Tolerance),  // Initialize tolerance with the specified value
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX)
{
    // Initialize randVect with createRandVect function
    createRandVector();

}


template <std::size_t N>
void MonteCarloSimulation<N>::simulatePhaseTransition() {
    int E = 0;
    int M = 0;
    float m = 0;
    float mExact = 1;
    float T = T_MIN;

    float error = 0;
    std::array<float, 2> prob;
    m = static_cast<float>(M) / N;
    while (T < T_MAX) {
        prob[0] = std::exp(-4 * lattice.getInteractionEnergy() / T);
        prob[1] = std::exp(-8 * lattice.getInteractionEnergy() / T);
        step = 0;
        mExact = mexact(T);
        error = std::abs(std::abs(m) - mExact);
        while (error > tolerance) {
            m = static_cast<float>(M) / N;
            error = std::abs(std::abs(m) - mExact);
            createRandVector();
            simulateStep(prob, lattice.getLattice(), random, M, E);
            step++;

        }

        T += Tstep;
        temperatures.pushback(T);
        monteCarloStepsResults.pushback(step);
        energyResults.pushback(E);
        magnetizationResults.pushback(m);
    }
}

template <std::size_t N>
void MonteCarloSimulation<N>::createRandVector() {
    for (int j = 0; j < N; j++) {
        randVect[j] = rand() % N;
    }
}


template <std::size_t N>
void MonteCarloSimulation<N>::translateMatrix(std::array<int, N>& inputMatrix) {
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
    // The memory managed by std::unique_ptr is automatically deallocated when it goes out of scope.
}

template <std::size_t N>
void MonteCarloSimulation<N>::flip(std::array<int, N>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) {
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

template <std::size_t N>
void MonteCarloSimulation<N>::simulateStep(std::array<float, 2> prob, std::array<int, N>& lattice, int& M, int& E) {
    for (unsigned long int i = 0; i < (N); i++) {
        int n = randVect[i];
        if (n != -1) {
            flip(lattice, prob, n, M, E);
        }
    }
}

template <std::size_t N>
void MonteCarloSimulation<N>::storeResultsToFile() const {
    // Open the file for writing
    std::ofstream outFile("result_" + std::to_string(N) + ".txt");

    // Check if the file is open
    if (!outFile.is_open()) {
        // Handle error: unable to open the file
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    // Write column headers
    outFile << "E M T N" << std::endl;

    // Determine the number of results to write
    std::size_t numResults = std::min({energyResults.size(), magnetizationResults.size(), temperatures.size(), monteCarloStepsResults.size()});

    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << energyResults[i] << " " << magnetizationResults[i] << " " << temperatures[i] << " " << monteCarloStepsResults[i] << std::endl;
    }

    // Close the file
    outFile.close();
}


