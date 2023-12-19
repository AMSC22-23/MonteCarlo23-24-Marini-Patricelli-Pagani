#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include "MontIntegration.hpp"
#include "Shape.hpp"
#include "HyperRectangle.hpp"
#include "HyperSphere.hpp"
#include <mpi.h>
#include "muParser.h"

int main(int argc, char** argv) {
    // Check the number of command line arguments
    if (argc != 4) {
        //@note: nice error, be more precise on the meaning
        std::cerr << "Usage: " << argv[0] << " <typeShape> <Nsamples> <paramsShapeFilename>" << std::endl;
        return -1;
    }

    MPI_Init(&argc, &argv);

    int r, s; //rank and size
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &s);

    // Parse command line arguments
    const auto typeShape = std::stoi(argv[1]);
    const auto N = std::strtol(argv[2], nullptr, 10);
    const std::string paramsShape = argv[3];

    // Create MontIntegration instance and shape
    MontIntegration montIntegration{};
    std::shared_ptr<Shape> shape;

    //@note: nice use of polymorphism
    switch (typeShape) {
        case 0: {
            shape = std::make_shared<HyperRectangle>(paramsShape);
            break;
        }
        case 1: {
            shape = std::make_shared<HyperSphere>(paramsShape);
            break;
        }
        default:
            std::cerr << "Invalid shape type." << std::endl;
            MPI_Finalize();
            return -1;
    }

    // Measure execution time
    auto start = std::chrono::system_clock::now();
    montIntegration.integrate(shape, N);
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    // Display results on rank 0
    if (r == 0) {
        std::cout << "Montecarlo Integral: " << montIntegration.getIntegral() << std::endl
                  << "Estimated error: " << std::sqrt(montIntegration.getVariance()) << std::endl
                  << "Target error: " << montIntegration.getIntegral() / std::sqrt(N) << std::endl
                  << "Shape dimension: " << shape->getDimensions() << std::endl
                  << "Shape volume: " << shape->getNorma() << std::endl
                  << "Elapsed time: " << elapsed_seconds.count() << " s" << std::endl;
    }

    MPI_Finalize();

    return 0;
}
