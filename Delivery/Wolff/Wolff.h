// MonteCarloSimulation.h

#ifndef MY_PROJECT_MONTE_CARLO_SIMULATION_H
#define MY_PROJECT_MONTE_CARLO_SIMULATION_H

#include <array>
#include <vector>
#include <memory>
#include <unordered_set>
#include <queue>
#include "SquareLattice.h"



class Wolff {
public:
    Wolff(float interactionStrength, int latticeSize,  float T_MIN, float T_MAX, float T_STEP , long int IT) ;
    void simulate_phase_transition() ;
    void store_results_to_file() const  ;
    



protected:
     void add_to_cluster(std::vector<int>& lattice, std::unordered_set<int>& cluster, std::queue<int>& spin_queue, float P, int i);
     void update(std::vector<int>& lattice, float T, const float J);
     float simulate(float T, std::vector<int>& lattice, float& energy, float& M, std::vector<int>& rand_vect, const float J);
     void createRandVect(std::vector<int>& rand_vect, int iterations);
     float calculate_magnetization_per_site(std::vector<int>& lattice);
     void create_rand_vect(std::vector<int>& rand_vect_0);


private:
    SquareLattice lattice;  // Use SquareLattice as a private member
    std::unique_ptr<std::vector<float> > MagnetizationResults; // Vector to store magnetization results.
    std::unique_ptr<std::vector<float> > Temperatures; // vector to store temperature visited
    float T_MIN;
    float T_MAX;
    float tolerance; // Tolerance value.
    float T_STEP;
    int L;
    int N ;
    long int IT;
};




#endif // MY_PROJECT_MONTE_CARLO_SIMULATION_H