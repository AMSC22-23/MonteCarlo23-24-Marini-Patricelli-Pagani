// MonteCarloSimulation.h

#ifndef MY_PROJECT_MONTE_CARLO_SIMULATION_H
#define MY_PROJECT_MONTE_CARLO_SIMULATION_H

#include <array>
#include <vector>
#include <memory>
#include "SquareLattice.h"
#include <mpi.h>



class Rem {
public:

    Rem(float interactionStrength, int latticeSize,  float T_MIN, float T_MAX, float T_STEP , long int IT) ;

void simulate_phase_transition(int argc, char* argv[]) ;
    void store_results_to_file() const  ;

protected:
    int find_set(int x, std::vector<int>& parent) ;
    void union_sets(int x, int y, std::vector<int>& parent, std::vector<int>& rank) ;
    void simulate_step(std::vector<int>& lattice, float P, std::vector<int>& parent, std::vector<int>& rank);
    void simulate(float P, std::vector<int>& lattice,int world_size, float& T, float J,int world_rank);
    float calculate_magnetization_per_site( std::vector<int>& lattice);
    void attempt_replica_exchange(std::vector<int>& lattice, float& T,float  J, int world_rank, int world_size, MPI_Comm comm);
    float evaluate(std::vector<int>& lattice, float J);

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
