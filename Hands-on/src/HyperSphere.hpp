#ifndef HYPERSPHERE_HPP
#define HYPERSPHERE_HPP

#include <mpi.h>

#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "Shape.hpp"

//@note: should put definition in source file
//       you were re-declaring the class in the source
//       which produces a conflict
//       Also you needed a semicolumn at the end

class HyperSphere : public Shape {
 public:
 //@note: should pass by const reference
  explicit HyperSphere(std::string filename) {
    std::ifstream inFile;
    inFile.open(filename, std::ios_base::in);

    double radius;

    inFile >> n;
    inFile >> function;
    inFile >> radius;

    center.reserve(n);
    center.resize(n);
    bounds.reserve(n);
    bounds.resize(n);

    if (radius == 0.0) exit(-1);

    for (int i = 0; i < n; ++i) {
      inFile >> center.at(i);
      bounds.at(i).x = center.at(i) - radius;
      bounds.at(i).y = center.at(i) + radius;
    }

    int r;  // rank
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    engine.seed(r);

    inFile.close();

    if (r == 0) calculateNorma();
  }
  std::vector<double> generateVector() override {
    std::vector<double> point;
    point.reserve(n);
    point.resize(n);

    sum = 0.0;

#pragma omp parallel for reduction(+ : sum)

    for (int j = 0; j < n; ++j) {
      std::uniform_real_distribution<double> distribution(bounds.at(j).x,
                                                          bounds.at(j).y);
      point.at(j) = distribution(engine);
      sum += (point.at(j) - center.at(j)) * (point.at(j) - center.at(j));
    }

    if (sum <= radius * radius) {
      return std::vector<double>();
    }

    return point;
  }

 private:
  double radius = 0.0;
  std::vector<double> center;
  std::vector<Edges> bounds;
  double sum;

  void calculateNorma() {
    double volTotale = 1.;

    volTotale *= std::pow(radius, n);
    volTotale *= std::pow(M_PI, (n / 2.));
    volTotale /= std::tgamma((n / 2.) + 1);

    norma = volTotale;
  }
};

#endif
