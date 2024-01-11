#include <chrono>
#include <vector>
#include <fstream> 
// Include the necessary header file
#include "Metropolis/DomainDecomposition/DomainDecomposition.h"
#include "Metropolis/SlidingWindow/SlidingWindow.h"
#include "Metropolis/CheckBoard/CheckBoard.h"
#include "Metropolis/SerialMetropolis/SerialMetropolis.h"


void store_performance_to_file(std::vector<float> time_results, std::vector<int> size, std::string filename){
    std::ofstream file;
    file.open("Performance_" + filename + ".txt");
    for(int i = 0; i < time_results.size(); i++){
        file << size[i] << " " << time_results[i] << std::endl;
    }
    file.close();
}

int main(){
    int L_MIN = 64;
    int L_MAX = 64;   
    float T_MIN=0.1;
    float T_MAX = 1.2;
    int T_STEP;
    float interactionStrength = 1;
    long int IT;
    std::vector<float> time_results;
    std::vector<int> size;
    int NUMTHREAD = 16;
    for(int L = L_MIN; L <= L_MAX; L *= 2){
        IT  = pow(L,4.4);
        std::cout<<"Simulation start for L = "<<L<<std::endl;
        CheckBoard simulation(interactionStrength, L, NUMTHREAD, T_MIN, T_MAX, T_STEP, IT);

        auto start = std::chrono::high_resolution_clock::now();

        simulation.simulate_phase_transition();

        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        time_results.emplace_back(duration.count());
        size.emplace_back(L*L);
        simulation.store_results_to_file();
        store_performance_to_file(time_results, size, "SlidingWindow");
    }
}


  

