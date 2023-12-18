#ifndef HYPERRECTANGLE_HPP
#define HYPERRECTANGLE_HPP

#include <mpi.h>

#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "Shape.hpp"

//@note: should put definition in source file
//       you were re-declaring the class in the source
//       which produces a conflict
//       Also you needed a semicolumn at the end

class HyperRectangle : public Shape {
 public:
  //@note: should pass by const reference
  explicit HyperRectangle(std::string filename) {
    std::ifstream inFile;
    inFile.open(filename, std::ios_base::in);

    inFile >> n;
    inFile >> function;
    bounds.reserve(n);
    bounds.resize(n);

    int i;

    for (i = 0; i < n; i++) {
      inFile >> bounds.at(i).x >> bounds.at(i).y;
      if (bounds.at(i).x == bounds.at(i).y) exit(-1);
      if (bounds.at(i).x > bounds.at(i).y)
        std::swap(bounds.at(i).x, bounds.at(i).y);
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    engine.seed(rank);

    inFile.close();

    if (rank == 0) calculateNorma();
  }

  std::vector<double> generateVector() override {
    std::vector<double> point;
    point.reserve(n);
    point.resize(n);

    for (int j = 0; j < n; ++j) {
      std::uniform_real_distribution<double> distribution(bounds.at(j).x,
                                                          bounds.at(j).y);
      point.at(j) = distribution(engine);
    }

    return point;
  }

 private:
  std::vector<Edges> bounds;
  std::string function;

  void calculateNorma() {
    double volTotale = norma;

    for (int i = 0; i < n; i++) {
      volTotale *= std::abs(bounds[i].x - bounds[i].y);
    }
    norma = volTotale;
  }
};

#endif
