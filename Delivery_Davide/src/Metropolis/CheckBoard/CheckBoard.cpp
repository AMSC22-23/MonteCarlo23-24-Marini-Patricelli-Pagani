
#include "CheckBoard.h"
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <random>


CheckBoard::CheckBoard(float interactionStrength, int latticeSize ,int NUMTHREAD , float T_MIN, float T_MAX, float T_STEP, long int IT)
: lattice(interactionStrength, latticeSize),  
        RandVect(),
        EnergyResults(),
        MagnetizationResults(),
        Temperatures(),
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        L(latticeSize),
        N(latticeSize * latticeSize),
        IT(IT),
        NUMTHREAD(NUMTHREAD),
        A(),
        ThreadStart(),
        NumFlipPerBlock(),
        rng(std::random_device{}()),
        dist(0.0, 1.0)
{   
    // Initialize tStart with std::make_unique
    ThreadStart = std::make_unique<std::vector<int> >();
    ThreadStart->reserve(NUMTHREAD);
    // Initialize A with set_block_size function
    A = set_block_size();
    NumFlipPerBlock = A*A;
    if (A == -1) {
    std::cerr << "Error: set_block_size failed\n";
    exit(EXIT_FAILURE);
    }
    omp_set_num_threads(NUMTHREAD);
    // Initialize RandVect with std::make_unique
    RandVect = std::make_unique<std::vector<int> >();
    RandVect->reserve(NumFlipPerBlock);
    // Initialize rand_vect with create_rand_vect function
    create_rand_vector();
    // Initialize EnergyResults with std::make_unique       
    EnergyResults = std::make_unique<std::vector<float> >();
    // Initialize MagnetizationResults with std::make_unique
    MagnetizationResults = std::make_unique<std::vector<float> >();
    // Initialize Temperatures with std::make_unique
    Temperatures = std::make_unique<std::vector<float> >();


}


void CheckBoard::simulate_phase_transition() {
    int E_loc ;
    int M_loc ;
    int deltaE ;
    int deltaM ;
    int step;
    float m = 0;
    float T = T_MIN;
    int M = lattice.getMagnetization();
    float error = 0;
    std::array<float, 2> prob;
    

#pragma omp parallel 
    {
        #pragma omp single nowait
        {
             while (T < T_MAX) {
                prob[0] = std::exp(-4 * lattice.getInteractionEnergy() / T);
                prob[1] = std::exp(-8 * lattice.getInteractionEnergy() / T);
                E_loc = 0;
                M_loc = 0; 
                deltaE = 0;
                deltaM = 0;
               for(int slide = 0; slide < ceil(IT/(NUMTHREAD*NumFlipPerBlock)); slide++ ){
                    create_rand_vector();
                    E_loc = 0;
                    M_loc = 0; 
                    //even blocks
                    for(int taskNum = 0; taskNum < NUMTHREAD; taskNum+=2){
                        #pragma omp task  firstprivate(M_loc, E_loc) shared(deltaM,deltaE)
                        {
                            simulate_step(prob,lattice.get_lattice(), M_loc, E_loc, (*ThreadStart)[taskNum]);
                            #pragma omp atomic update
                            deltaM += M_loc;
                            #pragma omp atomic update
                            deltaE += E_loc;
                        }
                    }
                    #pragma omp taskwait
                    //odd blocks
                    E_loc = 0;
                    M_loc = 0;
                    for(int taskNum = 1; taskNum < NUMTHREAD; taskNum+=2){
                        #pragma omp task  firstprivate(M_loc, E_loc) shared(deltaM,deltaE)
                        {
                            simulate_step(prob,lattice.get_lattice(), M_loc, E_loc, (*ThreadStart)[taskNum]);
                            #pragma omp atomic update
                            deltaM += M_loc;
                            #pragma omp atomic update
                            deltaE += E_loc;
                        }
                    }
                    #pragma omp taskwait
                }      
                  
                
            lattice.incrementMagnetization(deltaM);
            lattice.incrementEnergy(deltaE);
            m = static_cast<float>(lattice.getMagnetization()) / N;

            Temperatures->emplace_back(T);
            T += T_STEP;
            EnergyResults->emplace_back(1);
            MagnetizationResults->emplace_back(abs(m));
            lattice.restore_random_lattice();
             }
        }
    }
}


void CheckBoard::create_rand_vector() {
    //  N/NUMBLOCKS is the number of iteration for each task
    RandVect->clear();
    for (int j = 0;j<(NumFlipPerBlock);j++) {
        int r = static_cast<int>(dist(rng) * A);
        int c = static_cast<int>(dist(rng) * A);
        RandVect->emplace_back ((r * L + c));
        }
}

int CheckBoard::set_block_size() {
    int THREADPERSIDE = sqrt(NUMTHREAD);
    if(THREADPERSIDE*THREADPERSIDE != NUMTHREAD){
        std::cerr << "Error: Number of threads must be a perfect square." << std::endl;
        return -1;
    }
    else{
        if (L % THREADPERSIDE == 0) {
            int A_loc = L / THREADPERSIDE; // = larghezze di un blocco
        for(int i=0;i<  NUMTHREAD;i++){
                (*ThreadStart)[i] = (floor(i/THREADPERSIDE)*A_loc*L + i%THREADPERSIDE*A_loc); //thread at which each block start
            }
            return A_loc;
        }
        else {
            std::cerr << "Error: Unable to fill the line with the given number of blocks." << std::endl;
            return -1;
        }
    }
}



void CheckBoard::flip(std::vector<int>& lattice, std::array<float, 2>& prob, int site, int& M, int& E) {
    int sum = 0;

    if (site < L) {
        sum += lattice[site + L * (L - 1)];
    } else {
        sum += lattice[site - L];
    }
    if (site % L == 0) {
        sum += lattice[site + (L - 1)];
    } else {
        sum += lattice[site - 1];
    }

    if (site >= L * (L - 1)) {
        sum += lattice[site - L * (L - 1)];
    } else {
        sum += lattice[site + L];
    }
    if ((site + 1) % L == 0) {
        sum += lattice[site - (L - 1)];
    } else {
        sum += lattice[site + 1];
    }
    int delta = 2 * sum * lattice[site];
    if (delta <= 0) {
        lattice[site] = -lattice[site];
    } else if (delta == 4) {
        float rnd = dist(rng);
        if (rnd < prob[0]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    } else if (delta == 8) {
        float rnd = (dist(rng));
        if (rnd < prob[1]) {
            lattice[site] = -lattice[site];
        } else {
            return;
        }
    }
    M += 2 * lattice[site];
}


void CheckBoard::simulate_step (std::array<float, 2> prob, std::vector<int>& lattice, int& M_loc, int& E_loc, int offset) {
    int n;
    for (unsigned long int i = 0; i < (NumFlipPerBlock);i++) {
        n = (*RandVect)[i];
        flip(lattice, prob, n + offset, M_loc, E_loc);
    }
}


void CheckBoard::store_results_to_file() const {

    std::string filePath = "results/result_" + std::to_string(N) + ".txt";
    // Open the file for writing
    std::ofstream outFile("result_" + std::to_string(N) + ".txt");

    // Check if the file is open
    if (!outFile.is_open()) {
        // Handle error: unable to open the file
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    // Write column headers
    outFile << "E  M  T " << std::endl;

    // Determine the number of results to write
    std::size_t numResults = EnergyResults->size(); //they all have same lenght
    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << (*EnergyResults)[i] << " " << (*MagnetizationResults)[i] << " " << (*Temperatures)[i]  << std::endl;

    }

    // Close the file
    outFile.close();
}


