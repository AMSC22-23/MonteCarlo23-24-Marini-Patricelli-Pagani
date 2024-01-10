#include "MonteCarloSimulation.h"


int main(){
    int L ;
    int L_MAX = 4096;   
    float tollerance;
    for(L = 32; L <= L_MAX; L *= 2){
        tollerance = 2.0f / (L*L);
         std::cout<<"Simulation start for L = "<<L<<std::endl;
        MonteCarloSimulation simulation(1.0f, L, tollerance, 0.1f, 1.0f, 0.1f);
        simulation.simulatePhaseTransition();
        std::cout<<"Simulation done for L = "<<L<<std::endl;
        simulation.storeResultsToFile();
    }
    
    }

