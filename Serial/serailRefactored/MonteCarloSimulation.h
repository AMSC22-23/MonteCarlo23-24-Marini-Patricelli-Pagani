// MonteCarloSimulation.h

#ifndef MY_PROJECT_MONTE_CARLO_SIMULATION_H
#define MY_PROJECT_MONTE_CARLO_SIMULATION_H

#include <array>
#include <vector>
#include <memory>
#include "AbstractMonteCarloSimulation.h"
#include "SquareLattice.h"



class MonteCarloSimulation : public AbstractMonteCarloSimulation {
public:

    MonteCarloSimulation(float interactionStrength, int latticeSize,  float Tolerance, float T_MIN, float T_MAX, float T_STEP ) ;

    void simulatePhaseTransition() override;
    void storeResultsToFile() const override ;

protected:
    void createRandVector()  override
    {
         for (int j = 0; j < N; j++) {
        (*randVect)[j] = rand() % N;
         }
    }
    float mexact(float T) const override{
        return std::pow((1.0 - std::pow(std::sinh(2 * lattice.getInteractionEnergy() / T), -4)), 1.0 / 8.0);
    }
    void flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) override;
    void simulateStep(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E) override;

private:
    void translateMatrix(std::vector<int>& inputMatrix);
    SquareLattice lattice;  // Use SquareLattice as a private member
    std::unique_ptr<std::vector<int> > randVect; // Random vector of size N.
    std::unique_ptr<std::vector<float> > energyResults;   // Vector to store energy results.
    std::unique_ptr<std::vector<float> > magnetizationResults; // Vector to store magnetization results.
    std::unique_ptr<std::vector<int> > monteCarloStepsResults; //Vector to store number of monteCarlo step performed
    std::unique_ptr<std::vector<float> > temperatures; // vector to store temperature visited
    float T_MIN;
    float T_MAX;
    float tolerance; // Tolerance value.
    float T_STEP;
    const int L;
    const int N ;
};


#endif // MY_PROJECT_MONTE_CARLO_SIMULATION_H
