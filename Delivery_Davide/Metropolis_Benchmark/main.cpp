#include "Metropolis/SlidingWindow/SlidingWindow.h"
#include "Metropolis/DomainDecomposition/DomainDecomposition.h"
#include "Metropolis/CheckBoard/CheckBoard.h"
#include "Metropolis/Serial/SerialMetropolis.h"
#include <chrono>
#include <vector>
#include <fstream> // Include the necessary header file

void store_performance_to_file(std::vector<float> time_results, std::vector<int> size, std::string filename){
    std::ofstream file;
    file.open("Performance_" + filename + ".txt");
    for(int i = 0; i < time_results.size(); i++){
        file << size[i] << " " << time_results[i] << std::endl;
    }
    file.close();
}

int main(){
    int L ;
    int L_MAX = 128;   
    long int IT;
    std::vector<float> time_results;
    std::vector<int> size;
    int NUMTHREAD = 16;
    for(L = 64; L <= L_MAX; L *= 2){
        IT  = pow(L,4.4);
        std::cout<<"Simulation start for L = "<<L<<std::endl;
        SlidingWindow simulation(1.0f, L, NUMTHREAD, 0.1f, 1.2f, 0.1f, IT);

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


  

