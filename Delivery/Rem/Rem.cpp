#include "Rem.h"
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <random>
#include <omp.h> 
std::vector<std::mt19937> rngs;



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

//To correctly evaluate energy
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


float Rem::calculate_magnetization_per_site(std::vector<int>& lattice) {
    int total_spin = std::accumulate(lattice.begin(), lattice.end(), 0);
    return static_cast<float>(total_spin) / N;
}

int Rem::find_set(int x, std::vector<int>& parent) {
    while (x != parent[x]) {
        parent[x] = parent[parent[x]];  // Path compression
        x = parent[x];
    }
    return x;
}


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

