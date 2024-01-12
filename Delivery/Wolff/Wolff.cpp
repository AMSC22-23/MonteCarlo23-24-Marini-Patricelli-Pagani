#include "Wolff.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <random>
std::vector<std::mt19937> rngs;



Wolff::Wolff(float interactionStrength, int latticeSize,  float T_MIN, float T_MAX, float T_STEP, long int IT):
     lattice(interactionStrength, latticeSize),  
        MagnetizationResults(),
        Temperatures(),
        T_STEP(T_STEP) ,
        T_MIN(T_MIN),
        T_MAX(T_MAX),
        L(latticeSize),
        N(latticeSize * latticeSize),
        IT(IT)
{      

    // Initialize MagnetizationResults with std::make_unique
    MagnetizationResults = std::make_unique<std::vector<float> >();
    // Initialize temperatures with std::make_unique
    Temperatures = std::make_unique<std::vector<float> >();                     

}

void Wolff::simulate_phase_transition() {
    using namespace std;
    unsigned seed = time(0);
    srand(seed);
    float M = 0;
    float  energy =0;
    float T = T_MIN;
    using namespace std;
    std::vector<int> rand_vect;
    create_rand_vect(rand_vect);
    while (T < T_MAX) {
        simulate(T, lattice.get_lattice(),energy,  M,rand_vect,lattice.J);       
        Temperatures->push_back(T);
        T += T_STEP;
        M =  calculate_magnetization_per_site(lattice.get_lattice());
        MagnetizationResults->push_back(abs(M)); 
        lattice.restore_random_lattice();
        std::cout<<M<<std::endl;
    }

}

float Wolff::calculate_magnetization_per_site(std::vector<int>& lattice) {
    int total_spin = std::accumulate(lattice.begin(), lattice.end(), 0);
    return static_cast<float>(total_spin) / N;
}



void Wolff::add_to_cluster(std::vector<int>& lattice, std::unordered_set<int>& cluster, std::queue<int>& spin_queue, float P, int i) {
    if (cluster.find(i) == cluster.end()) {
        cluster.insert(i);
        std::vector<int> neighbors = {
            (i + L) % N,
            (i - L + N) % N,
            (i + 1) % L + (i / L) * L,
            (i - 1 + L) % L + (i / L) * L
        };
        for (int neighbor : neighbors) {
            if (lattice[neighbor] == lattice[i] && (rand() / (float)RAND_MAX) < P) {
                if (cluster.find(neighbor) == cluster.end()) {
                    spin_queue.push(neighbor);
                }
            }
        }
    }
}

void Wolff::update(std::vector<int>& lattice, float T,  float J) {
    float P = 1 - exp(-2 * J / T);
    int seed = rand() % N;
    std::unordered_set<int> cluster;
    std::queue<int> spin_queue;
    spin_queue.push(seed);
    while (!spin_queue.empty()) {
        int current = spin_queue.front();
        spin_queue.pop();
        add_to_cluster(lattice, cluster, spin_queue, P, current);
    }
    for (int i : cluster) {
        lattice[i] *= -1;
    }
}


float Wolff::simulate(float T,std::vector <int> & lattice, float& energy, float& M,std::vector<int>& rand_vect,  float J) {
    using namespace std;
    vector<float> prob(2);
    prob[0] = exp(-4 * J/T);
    prob[1] = exp(-8 * J/T);
    vector<float> energy_vec(1);
    energy_vec[0] =  energy;
    vector<float> m(1);
    m[0] = (float)M/N;
    vector<int> t_axis(1);
    t_axis[0] = 0;
    int n;
    for (int i = 1; i < IT; i++) {
        update(lattice, T, J);
        if (i % N == 0) {
            m.push_back((float) M / N);
            energy_vec.push_back(energy);
            t_axis.push_back(i / N);
        }
    }
    return m.back();
}

void Wolff::create_rand_vect(std::vector<int>& rand_vect_0) {
    int i;
    for (i = 0; i < IT; i++) {
        rand_vect_0.push_back(rand() % N );
    }
}




void Wolff::store_results_to_file() const {
    // Open the file for writing
    std::ofstream outFile("result_" + std::to_string(N) + ".txt");

    // Check if the file is open
    if (!outFile.is_open()) {
        // Handle error: unable to open the file
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
        return;
    }

    // Write column headers
    outFile << " M  T" << std::endl;

    // Determine the number of results to write
    std::size_t numResults = MagnetizationResults->size(); //they all have same lenght
    // Write results to the file
    for (std::size_t i = 0; i < numResults; ++i) {
        // Write data for each row
        outFile << ( (*MagnetizationResults))[i] << " " << (*Temperatures)[i] << std::endl;

    }

    // Close the file
    outFile.close();
}


