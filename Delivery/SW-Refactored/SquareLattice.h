#ifndef SQUARELATTICE_H
#define SQUARELATTICE_H

#include "AbstractLattice.h"
#include <iostream>
#include <memory>
#include <vector>

class SquareLattice : public AbstractLattice {
public:
    SquareLattice(float interactionStrength, int latticeSize);

    void print_lattice() const ;
    float evaluate_energy() const ;
    void initialize() override;
    std::vector<int>& get_lattice();
    float get_interaction_energy() const override {
        return J;
    }

    void restore_random_lattice() {
        lattice->assign(randomLattice->begin(), randomLattice->end());
        M = M_rand;
        E = E_rand;
    }

private:
    const int L;
    const int J;
    int N; 
    int M;
    float E;
    int M_rand;
    float E_rand;
    std::unique_ptr<std::vector<int>> lattice;
    std::unique_ptr<std::vector<int>> randomLattice;
};


#endif // SQUARELATTICE_H
