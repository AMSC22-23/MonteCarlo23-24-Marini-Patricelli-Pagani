#ifndef ABSTRACT_MONTE_CARLO_SIMULATION_H
#define ABSTRACT_MONTE_CARLO_SIMULATION_H

#include <array>

class AbstractMonteCarloSimulation {
public:
    virtual void simulatePhaseTransition() = 0;
protected:
    virtual void createRandVector() = 0;
    virtual void flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) = 0;
    virtual void simulateStep(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E) = 0;
    virtual void storeResultsToFile() const = 0;
    virtual float mexact(float T) const = 0;

    
};

#endif // ABSTRACT_MONTE_CARLO_SIMULATION_H
