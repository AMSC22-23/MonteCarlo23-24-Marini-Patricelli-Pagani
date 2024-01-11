#ifndef SQUARELATTICE_H
#define SQUARELATTICE_H

#include "AbstractLattice.h"
#include <iostream>
#include <memory>
#include <vector>

class SquareLattice : public AbstractLattice {
public:
    SquareLattice(float interactionStrength, int latticeSize);

    void printLattice() const override;
    float evaluateEnergy() const override;
    void initialize() override;
    std::vector<int>& get_lattice();
    float getInteractionEnergy() const override {
        return J;
    }
    float getMagnetization() const {
        return M;
    }
    void incrementMagnetization(int Madd) {
        M += Madd;
    }
    void incrementEnergy(int Eadd) {
        E += Eadd;
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
