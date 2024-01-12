#ifndef ABSTRACT_MONTE_CARLO_SIMULATION_H
#define ABSTRACT_MONTE_CARLO_SIMULATION_H

#include <array>

class AbstractMonteCarloSimulation {
public:
    virtual void simulate_phase_transition() = 0;
protected:
    virtual void create_rand_vector() = 0;
    virtual void simulate_step(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E) = 0;
    virtual void store_results_to_file() const = 0;
    virtual float mexact(float T) const = 0;

    
};

#endif // ABSTRACT_MONTE_CARLO_SIMULATION_H
