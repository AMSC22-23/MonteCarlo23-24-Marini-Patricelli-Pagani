#include "Rem.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <random>
#include <omp.h> 
std::vector<std::mt19937> rngs;

/**
 * @brief Constructor for the Rem class.
 * 
 * Initializes various parameters and sets up data structures for the simulation.
 * 
 * @param interactionStrength The strength of the interaction in the lattice.
 * @param latticeSize The size of the lattice.
 * @param T_MIN Minimum temperature for the simulation.
 * @param T_MAX Maximum temperature for the simulation.
 * @param T_STEP Step size for temperature increment.
 * @param IT Number of iterations for the simulation.
 */

Rem::Rem(float interactionStrength, int latticeSize,  float T_MIN, float T_MAX, float T_STEP, long int IT):
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

/**
 * @brief Simulates the phase transition of the lattice.
 * 
 * This function orchestrates the overall simulation, managing temperature changes,
 * lattice updates, and data collection.
 * 
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * */

void Rem::simulate_phase_transition(int argc, char* argv[]) {
    
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand(static_cast<unsigned>(time(nullptr)) + world_rank);
 
    float m = 0;
    float T = T_MIN;
    using namespace std;
    // Set up the number of threads for OpenMP
    omp_set_num_threads(omp_get_max_threads());
    int num_threads = omp_get_max_threads();
    rngs.resize(num_threads); // Seed each RNG with a unique seed

    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    for (int i = 0; i < num_threads; ++i) {
        rngs[i].seed(seed + i);
    }
    
    while (T < T_MAX) {
        float P = 1 - exp(-2 * lattice.get_interaction_energy() / T);
        simulate(P, lattice.get_lattice(),world_size,T,lattice.J,world_rank);       
        Temperatures->push_back(T);
        T += T_STEP;
        m =  calculate_magnetization_per_site(lattice.get_lattice());
        MagnetizationResults->push_back(abs(m)); 
        lattice.restore_random_lattice();
        std::cout<<m<<std::endl;
    }
}

/**
 * @brief Evaluates, using OpenMP parallelism, the energyof the given lattice configuration.
 * 
 * This function calculates the total energy of the lattice based on its current state
 * and the interaction strength J.
 * 
 * @param lattice The vector representing the lattice configuration.
 * @param J The interaction strength.
 * @return float The calculated energy of the lattice.
 */

float Rem::evaluate(std::vector<int>& lattice, float J) {
    int sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i++) {
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L] * lattice[i] * 2;
        }
        if (i % L != 0) { //NO FIRST COLUMN
            sum += lattice[i - 1] * lattice[i] * 2;
        }

        if (i >= L * (L - 1)) { //LAST ROW
            sum += lattice[i - L * (L - 1)] * lattice[i] * 2; //times 2 two count also contribute where i = 0
        }
        if ((i + 1) % L == 0) { //LAST COLUMN
            sum += lattice[i - (L - 1)] * lattice[i] * 2;
        }
    }
    return -J * sum;
}

/**
 * @brief Calculates the magnetization per site of the lattice.
 * 
 * This function computes the average magnetization across the lattice.
 * 
 * @param lattice The vector representing the lattice configuration.
 * @return float The average magnetization per site.
 */


float Rem::calculate_magnetization_per_site(std::vector<int>& lattice) {
    int total_spin = std::accumulate(lattice.begin(), lattice.end(), 0);
    return static_cast<float>(total_spin) / N;
}

/**
 * @brief Finds the set of the given element in the disjoint-set data structure.
 * 
 * Implements the "find" operation of the union-find algorithm with path compression.
 * 
 * @param x The element whose set is to be found.
 * @param parent The vector representing the parent of each element.
 * @return int The root of the set containing 'x'.
 */

int Rem::find_set(int x, std::vector<int>& parent) {
    while (x != parent[x]) {
        parent[x] = parent[parent[x]];  // Path compression
        x = parent[x];
    }
    return x;
}

/**
 * @brief Merges the sets containing the two given elements.
 * 
 * Implements the "union" operation of the union-find algorithm.
 * 
 * @param x The first element.
 * @param y The second element.
 * @param parent The vector representing the parent of each element.
 * @param rank The vector representing the rank of each element.
 */

void Rem::union_sets(int x, int y, std::vector<int>& parent, std::vector<int>& rank) {
    x = find_set(x, parent);
    y = find_set(y, parent);
    if (x != y) {
        if (rank[x] < rank[y]) {
            std::swap(x, y);
        }
        parent[y] = x;
        if (rank[x] == rank[y]) {
            rank[x]++;
        }
    }
}

/**
 * @brief Performs a single simulation step.
 * 
 * This function updates the lattice configuration for a single step of the simulation.
 * 
 * @param lattice The vector representing the lattice configuration.
 * @param P The probability of forming a bond between aligned spins.
 * @param parent The vector representing the parent of each element in the disjoint set.
 * @param rank The vector representing the rank of each element in the disjoint set.
 */

void Rem::simulate_step(std::vector<int>& lattice, float P, std::vector<int>& parent, std::vector<int>& rank) {
    

    // Reset parent and rank for each site
    std::fill(parent.begin(), parent.end(), 0);
    std::fill(rank.begin(), rank.end(), 0);
    std::iota(parent.begin(), parent.end(), 0); // Set each site as its own parent initially

    // Form clusters
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int thread_id = omp_get_thread_num();
        std::mt19937& local_rng = rngs[thread_id];
        std::uniform_real_distribution<float> dist(0.0, 1.0);

        int right = (i % L == L - 1) ? i - (L - 1) : i + 1;
        int down = (i + L) % N;

        if (lattice[i] == lattice[right] && dist(local_rng) < P) {
            union_sets(i, right, parent, rank);
        }
        if (lattice[i] == lattice[down] && dist(local_rng) < P) {
            union_sets(i, down, parent, rank);
        }
    }

    std::vector<int> flip_decision(N);
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int thread_id = omp_get_thread_num();
        std::mt19937& local_rng = rngs[thread_id];
        std::uniform_int_distribution<int> dist(0, 1);

        flip_decision[i] = dist(local_rng);
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        if (flip_decision[find_set(i, parent)]) {
            lattice[i] *= -1;
        }
    }
}

/**
 * @brief Conducts the simulation over multiple iterations.
 * 
 * This function manages the simulation process, including updating lattice states
 * and attempting replica exchanges.
 * 
 * @param P The probability of forming a bond between aligned spins.
 * @param lattice The vector representing the lattice configuration.
 * @param world_size The size of the MPI world (number of processes).
 * @param T The current temperature.
 * @param J The interaction strength.
 * @param world_rank The rank of the current MPI process.
 */


void Rem::simulate(float P, std::vector<int>& lattice,int world_size, float& T, float J,int world_rank) {
    using namespace std;
    vector<int> parent(N), rank(N, 0);
    for (int i = 1; i < IT; ++i) {
        simulate_step(lattice, P, parent, rank);
                 // Attempt to exchange configuration with neighboring replica every 5 iterations
        if (i % 5 == 0 && world_size > 1) {
            attempt_replica_exchange(lattice, T, J,world_rank, world_size, MPI_COMM_WORLD);
        }
    }

}

/**
 * @brief Attempts to exchange lattice configurations with a neighboring replica.
 * 
 * This function implements the replica exchange mechanism in parallel simulations.
 * 
 * @param lattice The vector representing the lattice configuration.
 * @param T The current temperature of the replica.
 * @param J The interaction strength.
 * @param world_rank The rank of the current MPI process.
 * @param world_size The size of the MPI world.
 * @param comm The MPI communicator.
 */



void Rem::attempt_replica_exchange(std::vector<int>& lattice, float& T, float J, int world_rank, int world_size, MPI_Comm comm) {
    int partner;
    if (world_rank % 2 == 0) {
        partner = world_rank + 1;
    } else {
        partner = world_rank - 1;
    }

    // Boundary condition for the replicas
    if (partner < 0 || partner >= world_size) return;

    float my_energy = evaluate(lattice,J);
    float partner_energy;
    float my_T = T, partner_T;

    // Exchange temperatures with the partner
    MPI_Sendrecv(&my_T, 1, MPI_FLOAT, partner, 0, &partner_T, 1, MPI_FLOAT, partner, 0, comm, MPI_STATUS_IGNORE);

    // Only exchange energies and attempt swapping if temperatures are different
    if (my_T != partner_T) {
        MPI_Sendrecv(&my_energy, 1, MPI_FLOAT, partner, 1, &partner_energy, 1, MPI_FLOAT, partner, 1, comm, MPI_STATUS_IGNORE);

        // Metropolis criterion
        float delta = (1.0/my_T - 1.0/partner_T) * (partner_energy - my_energy);
        if (delta < 0 || exp(-delta) > static_cast<float>(rand()) / RAND_MAX) {
            // Accept the swap
            std::vector<int> partner_lattice(N);
            MPI_Sendrecv(lattice.data(), N, MPI_INT, partner, 2, partner_lattice.data(), N, MPI_INT, partner, 2, comm, MPI_STATUS_IGNORE);
            lattice.swap(partner_lattice);
        }
    }
}


/**
 * @brief Stores the results of the simulation to a file.
 * 
 * This function writes the magnetization results and corresponding temperatures
 * to a file for later analysis.
 */

//With MPI parallelism, I am not sure this always store the correct result, the result corresponding o the requested replica.

void Rem::store_results_to_file() const {
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

