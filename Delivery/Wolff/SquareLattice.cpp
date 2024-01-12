#include "SquareLattice.h"

SquareLattice::SquareLattice(float interactionStrength, int latticeSize) 
    : AbstractLattice(),
     L(latticeSize),
     N(latticeSize * latticeSize),
     J(interactionStrength), 
     lattice(),
     randomLattice(),
     M_rand(),
     E_rand(),
     M(),
     E() 
{
    lattice = std::make_unique<std::vector<int>>(N);
    randomLattice = std::make_unique<std::vector<int>>(N);
    initialize();
    restore_random_lattice();
}

void SquareLattice::print_lattice() const {
    for (int i = 0; i < N; i++) {
        if (i % L == 0)
            std::cout << std::endl;
        if ((*lattice)[i] == -1)
            std::cout << "o ";
        else
            std::cout << "x ";
    }
    std::cout << std::endl;
}


float SquareLattice::evaluate_energy() const {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        if (i >= L) // NO FIRST ROW
            sum += (*lattice)[i - L] * (*lattice)[i] * 2;
        if (i % L != 0) // NO FIRST COLUMN
            sum += (*lattice)[i - 1] * (*lattice)[i] * 2;
        if (i >= L * (L - 1)) // LAST ROW
            sum += (*lattice)[i - L * (L - 1)] * (*lattice)[i] * 2;
        if ((i + 1) % L == 0) // LAST COLUMN
            sum += (*lattice)[i - (L - 1)] * (*lattice)[i] * 2;
    }
    return -J * sum;
}


void SquareLattice::initialize() {
    int sum = 0;
    M_rand = 0;

    for (int i = 0; i < N; i++) {
        int k = rand() % 2;
        if (k == 0) {
            (*randomLattice)[i] = -1;
            M_rand -= 1;
        } else {
            (*randomLattice)[i] = 1;
            M_rand += 1;
        }

        if (i >= L) // NO FIRST ROW
            sum += (*randomLattice)[i - L] * (*randomLattice)[i] * 2;
        if (i % L != 0) // NO FIRST COLUMN
            sum += (*randomLattice)[i - 1] * (*randomLattice)[i] * 2;
        if (i >= L * (L - 1)) // LAST ROW
            sum += (*randomLattice)[i - L * (L - 1)] * (*randomLattice)[i] * 2;
        if ((i + 1) % L == 0) // LAST COLUMN
            sum += (*randomLattice)[i - (L - 1)] * (*randomLattice)[i] * 2;
    }

    E_rand = -J * sum;
}



std::vector<int>& SquareLattice::get_lattice() {
    return *lattice;
}


