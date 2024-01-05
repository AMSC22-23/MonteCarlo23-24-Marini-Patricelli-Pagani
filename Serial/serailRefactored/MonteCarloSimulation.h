// MonteCarloSimulation.h

#ifndef MY_PROJECT_MONTE_CARLO_SIMULATION_H
#define MY_PROJECT_MONTE_CARLO_SIMULATION_H

#include <array>
#include <vector>
#include <memory>
#include "AbstractMonteCarloSimulation.h"
#include "SquareLattice.h"



template <std::size_t N>
class MonteCarloSimulation : public AbstractMonteCarloSimulation<N> {
public:

    MonteCarloSimulation(float interactionStrength, int latticeSize,  float Tolerance, float T_MIN, float T_MAX, float T_STEP ) ;

    void simulatePhaseTransition();
    void storeResultsToFile() const

protected:
    void createRandVector() override;
    void flip(std::array<int, N>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) override;
    void simulateStep(std::array<float, 2> prob, std::array<int, N>& lattice, int& M, int& E) override;

private:
    void translateMatrix(std::array<int, N>& inputMatrix);
    SquareLattice<N> lattice;  // Use SquareLattice as a private member
    std::unique_ptr<std::array<int, N>> randVect; // Random vector of size N.
    std::unique_ptr<std::vector<float>> energyResults;   // Vector to store energy results.
    std::unique_ptr<std::vector<float>> magnetizationResults; // Vector to store magnetization results.
    std::unique_ptr<std::vector<int>> monteCarloStepsResults; //Vector to store number of monteCarlo step performed
    std::unique_ptr<std::vector<float>> temperatures; // vector to store temperature visited
    float T_MIN;
    float T_MAX;
    float tolerance; // Tolerance value.
    float T_STEP;

};



#include "MonteCarloSimulation.tpp"  // Include the source file

#endif // MY_PROJECT_MONTE_CARLO_SIMULATION_H
