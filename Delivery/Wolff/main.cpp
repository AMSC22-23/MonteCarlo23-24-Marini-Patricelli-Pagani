#include "Wolff.h"
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <unordered_set>
#include <queue>
#include <map>
#include <numeric>
#include <omp.h>
#include <random>


int main() {
    int L_MIN;
    int L_MAX;
    float T_MIN;
    float T_MAX;
    float T_STEP;
    float J;

    std::cout << "Enter L_MIN: ";
    std::cin >> L_MIN;

    std::cout << "Enter L_MAX: ";
    std::cin >> L_MAX;

    std::cout << "Enter T_MIN: ";
    std::cin >> T_MIN;

    std::cout << "Enter T_MAX: ";
    std::cin >> T_MAX;

    std::cout << "Enter T_STEP: ";
    std::cin >> T_STEP;

    std::cout << "Enter J: ";
    std::cin >> J;

    std::vector<float> ITVector;

    float IT;
    for (int L = L_MIN; L <= L_MAX; L *= 2) {
        std::cout << "Enter IT for L = " << L << ": ";
        std::cin >> IT;
        ITVector.emplace_back(IT);
    }

    int L = L_MIN;
    for(int i = 0; i < ITVector.size(); i++) {
         std::cout<<"Simulation start for L = "<<L<<std::endl;
        Wolff simulation(J, L, T_MIN, T_MAX, T_STEP, ITVector[i]);
        simulation.simulate_phase_transition();
        std::cout<<"Simulation done for L = "<<L<<std::endl;
        simulation.store_results_to_file();
        L *= 2;
    }
    
    }
