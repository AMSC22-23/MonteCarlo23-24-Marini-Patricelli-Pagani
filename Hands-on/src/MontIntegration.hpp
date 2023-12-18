#pragma once

//@note: here you are declaring a class, must be a header file

#include <functional>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <mpi.h>
#include "Shape.hpp"
#include "muParser.h"

class MontIntegration {
public:
    MontIntegration() = default;

private:
    double integral;
    double variance;

//@note: this needs to be public
public:
    void integrate(std::shared_ptr<Shape> shape, long N) {
        int r,s;  //rank e size

        MPI_Comm_rank(MPI_COMM_WORLD, &r);
        MPI_Comm_size(MPI_COMM_WORLD, &s);

        MPI_Bcast(&N, 1, MPI_LONG, 0, MPI_COMM_WORLD);

        const long points = N / s + (r < (N % s));
        std::vector<double> point(shape->getDimensions());
        std::vector<double> sample;

        // Parser initialization
        mu::Parser p;
        p.SetExpr(shape->getFunction());

        double sum = 0.0;
        double sum_2= 0.0;

        //@note: never hard code the number of threads
        #pragma omp parallel for num_threads(12) default(none) \
            shared(shape, N, r, p, points) \
            reduction(+ : sum, sum_2) private(point, sample)
        for (int i = 0; i < points; ++i) {
            sample = shape->generateVector();

            for (int j = 0; j < shape->getDimensions(); ++j) {
                std::string num_0 = "x";
                std::string num = std::to_string(j);
                num_0 = num_0 + num;
                p.DefineVar(num_0, &sample.at(j));
            }

            auto fi = p.Eval();

            sum += fi;
            sum_2 += fi * fi;
        }

        double sumTot = 0.0; 
        double sumTot_2 = 0.0;
        MPI_Reduce(&sum, &sumTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&sum_2, &sumTot_2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (r == 0) {
            integral = shape->getNorma() * sumTot / N;
            variance = shape->getNorma() * shape->getNorma() * ((sumTot_2 - (sumTot * sumTot) / N) / (N - 1)) / N;

            }
    };

    double getIntegral() const { return integral; }
    double getVariance() const { return variance; }    
};
