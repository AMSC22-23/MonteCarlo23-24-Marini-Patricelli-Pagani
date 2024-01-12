/*The structure of this CUDA program is specifically designed for development on cloud-based platforms like Google Colab. 
These platforms were chosen for their suitability in handling CUDA programming, especially for those without access to local CUDA hardware. 
The streamlined code structure, avoiding traditional `.cuh` and `.cu` files and `CMake`, is more practical for these environments,
 focusing on ease of use and efficiency in building and executing the program.
*/

#include <cuda.h>
#include <curand_kernel.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

/* 

This CUDA program is designed with flexibility in mind, allowing for key parameters to be adjusted according to the specific requirements of the simulation.
 The primary configurable parameters include:

L (Lattice Size): Determines the dimensions of the lattice (N = L * L). It's important that L is a multiple of two for proper lattice structure.

J (Interaction Strength): This parameter influences the interaction strength within the system, and can be modified to study various physical behaviors.

NTHREADS (Number of GPU Threads): Critical for performance optimization. Choose a value that allows the GPU to efficiently manage computation.

IT (Number of Iterations): Defines the total number of iterations for the simulation, crucial for ensuring the accuracy and convergence of the results.

Performance and Convergence Guidelines

Thread Count and Lattice Size: It's recommended that L be divisible by NTHREADS to ensure efficient workload distribution and to avoid indexing issues. 
If you modify L, adjust NTHREADS accordingly to maintain this divisibility.

Optimal Settings: The program was tested with L = 1024, IT = 2e9, and NTHREADS = 256. These settings are known to provide a high confidence of convergence in simulations.
Adjusting for Larger Lattices: If you choose to increase L, be mindful that it might require a proportional adjustment in NTHREADS and IT depending on your GPU capabilities. 
A larger L generally necessitates more iterations to ensure convergence, and the thread count may need to be modified to maintain optimal performance.

Convergence Assurance: With L = 1024 and IT = 2e9, the program is highly likely to achieve reliable convergence, making these settings a good starting point for simulations. 

*/

#define L 1024
#define N (L*L)
#define J 1.00
#define IT 5e8 // Number of iterations, should be divisible by 2 for even updates
#define NTHREADS 256 // Number of GPU threads

__device__ int get_index(int row, int col);
__device__ int delta_energy(bool* lattice, int r, int c);
__global__ void flip_spins(bool* lattice, float* prob, float* energy, curandState* states, bool update_black);
__global__ void setup_rand_kernel(curandState* state, unsigned long seed);
__global__ void initialize_lattice_kernel(bool* lattice, curandState* states);

int main() {
    bool* dev_lattice;
    cudaMalloc((void**)&dev_lattice, N * sizeof(bool));

    curandState* dev_states;
    cudaMalloc((void**)&dev_states, N * sizeof(curandState));

    dim3 blocksPerGrid((N + NTHREADS - 1) / NTHREADS, 1, 1);
    dim3 threadsPerBlock(NTHREADS, 1, 1);

    unsigned long seed = static_cast<unsigned long>(time(nullptr));
    setup_rand_kernel << <blocksPerGrid, threadsPerBlock >> > (dev_states, seed);

    float* dev_energy;
    cudaMalloc((void**)&dev_energy, sizeof(float));

    float energy = 0.0f;
    cudaMemcpy(dev_energy, &energy, sizeof(float), cudaMemcpyHostToDevice);

    initialize_lattice_kernel << <blocksPerGrid, threadsPerBlock >> > (dev_lattice, dev_states);

    float* dev_probabilities;
    cudaMalloc((void**)&dev_probabilities, 2 * sizeof(float));

    for (float T = 0.2f; T <= 3.0f; T += 0.1f) {
        clock_t start_time = clock();

        float prob[2] = { exp(-4 * J / T), exp(-8 * J / T) };
        cudaMemcpy(dev_probabilities, prob, 2 * sizeof(float), cudaMemcpyHostToDevice);

        for (unsigned long i = 0; i < IT / N; i += 2) {
            flip_spins << <blocksPerGrid, threadsPerBlock >> > (dev_lattice, dev_probabilities, dev_energy, dev_states, true);
            flip_spins << <blocksPerGrid, threadsPerBlock >> > (dev_lattice, dev_probabilities, dev_energy, dev_states, false);
        }

        cudaDeviceSynchronize();

        clock_t end_time = clock();

        double elapsed_secs = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;

        cudaMemcpy(&energy, dev_energy, sizeof(float), cudaMemcpyDeviceToHost);

        std::cout << "Temperature: " << T << std::endl;
        std::cout << "Final Energy: " << energy / N << std::endl;
        std::cout << "Simulation time (seconds): " << elapsed_secs << std::endl << std::endl;
    }

    cudaFree(dev_lattice);
    cudaFree(dev_states);
    cudaFree(dev_energy);
    cudaFree(dev_probabilities);

    return 0;
}



__device__ int get_index(int row, int col) {
    return (row * L + col) % N;
}

// The lattice uses boolean values, true for spin up (equivalent to 1) and false for spin down (equivalent to -1)
__device__ int delta_energy(bool* lattice, int r, int c) {
    int sum = lattice[get_index((r - 1 + L) % L, c)]
        + lattice[get_index((r + 1) % L, c)]
        + lattice[get_index(r, (c - 1 + L) % L)]
        + lattice[get_index(r, (c + 1) % L)];
    sum = 2 * sum - 4; // Convert sum from [0, 4] to [-4, 4] to match the original spin values
    int spin = lattice[get_index(r, c)] ? 1 : -1; // Convert bool to equivalent spin value
    return 2 * spin * sum;
}

__global__ void flip_spins(bool* lattice, float* prob, float* energy, curandState* states, bool update_black) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    int r = idx / L;
    int c = idx % L;
    bool is_black = ((r + c) % 2 == 0);

    if (is_black == update_black) {
        int delta = delta_energy(lattice, r, c);
        float rnd = curand_uniform(&states[idx]);

        if (delta <= 0 || (delta == 4 && rnd < prob[0]) || (delta == 8 && rnd < prob[1])) {
            lattice[get_index(r, c)] = !lattice[get_index(r, c)];
            atomicAdd(energy, delta * J);
            // Removed magnetization update
        }
    }
}

__global__ void initialize_lattice_kernel(bool* lattice, curandState* states) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        float randVal = curand_uniform(&states[idx]);
        lattice[idx] = (randVal < 0.5f);
        // Removed magnetization calculation and update
    }
}

__global__ void setup_rand_kernel(curandState* state, unsigned long seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        curand_init(seed, idx, 0, &state[idx]);
    }
}
