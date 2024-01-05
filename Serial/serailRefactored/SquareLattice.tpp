#ifndef SQUARELATTICE_TPP
#define SQUARELATTICE_TPP

template <int N>
SquareLattice<N>::SquareLattice(float interactionStrength, int latticeSize) 
    : AbstractLattice<N>(), L(latticeSize), lattice()
{
    initialize();
}

template <int N>
void SquareLattice<N>::printLattice() const {
    for (int i = 0; i < N; i++) {
        if (i % L == 0)
            std::cout << std::endl;
        if (lattice[i] == -1)
            std::cout << "o ";
        else
            std::cout << "x ";
    }
    std::cout << std::endl;
}

template <int N>
float SquareLattice<N>::evaluateEnergy() const {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        if (i >= L) // NO FIRST ROW
            sum += lattice[i - L] * lattice[i] * 2;
        if (i % L != 0) // NO FIRST COLUMN
            sum += lattice[i - 1] * lattice[i] * 2;
        if (i >= L * (L - 1)) // LAST ROW
            sum += lattice[i - L * (L - 1)] * lattice[i] * 2;
        if ((i + 1) % L == 0) // LAST COLUMN
            sum += lattice[i - (L - 1)] * lattice[i] * 2;
    }
    return -J * sum;
}

template <int N>
void SquareLattice<N>::initialize() {
    int sum = 0;
    int M = 0;

    for (int i = 0; i < N; i++) {
        int k = rand() % 2;
        if (k == 0) {
            lattice[i] = -1;
            M -= 1;
        } else {
            lattice[i] = 1;
            M += 1;
        }

        if (i >= L) // NO FIRST ROW
            sum += lattice[i - L] * lattice[i] * 2;
        if (i % L != 0) // NO FIRST COLUMN
            sum += lattice[i - 1] * lattice[i] * 2;
        if (i >= L * (L - 1)) // LAST ROW
            sum += lattice[i - L * (L - 1)] * lattice[i] * 2;
        if ((i + 1) % L == 0) // LAST COLUMN
            sum += lattice[i - (L - 1)] * lattice[i] * 2;
    }

    evaluateEnergy();
}


template <int N>
const std::array<int, N>& SquareLattice<N>::getLattice() const {
    return lattice;
}


#endif // SQUARELATTICE_TPP
