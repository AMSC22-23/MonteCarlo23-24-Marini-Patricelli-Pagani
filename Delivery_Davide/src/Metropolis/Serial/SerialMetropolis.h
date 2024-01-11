// SerialMetropolis.h

#ifndef MY_PROJECT_MONTE_CARLO_SIMULATION_H
#define MY_PROJECT_MONTE_CARLO_SIMULATION_H

#include <array>
#include <vector>
#include <memory>
#include "../AbstractMonteCarloSimulation.h"
#include "../../Lattice/SquareLattice.h"



class SerialMetropolis : public AbstractMonteCarloSimulation {
public:

    SerialMetropolis(float interactionStrength, int latticeSize , float T_MIN, float T_MAX, float T_STEP , long int IT) ;

    void simulate_phase_transition() override;
    void store_results_to_file() const override ;

protected:
    void create_rand_vector()  override;
    void flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E);
    void simulate_step(std::array<float, 2> prob, std::vector<int>& lattice, int& M, int& E, int offset = 0) override;

private:
    SquareLattice lattice;  // Use SquareLattice as a private member
    std::unique_ptr<std::vector<int> > RandVect; // Random vector of size N.
    std::unique_ptr<std::vector<float> > EnergyResults;   // Vector to store energy results.
    std::unique_ptr<std::vector<float> > MagnetizationResults; // Vector to store magnetization results. 
    std::unique_ptr<std::vector<float> > Temperatures; // vector to store temperature visited
    std::unique_ptr<std::vector<int> > ThreadStart; // vector to store monte carlo steps
    float T_MIN;
    float T_MAX;
    float T_STEP;
    const int L;
    const int N ;
    const long int IT;
};


#endif // MY_PROJECT_MONTE_CARLO_SIMULATION_H
