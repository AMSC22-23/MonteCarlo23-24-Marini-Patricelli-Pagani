// AbstractLattice.h

#ifndef ABSTRACTLATTICE_H
#define ABSTRACTLATTICE_H

/**
 * @brief Abstract base class for a lattice structure.
 * 
 * This class defines an interface for various lattice implementations.
 * It includes methods for initialization, energy evaluation, printing
 * the lattice, and getting interaction energy.
 */

class AbstractLattice {
public:

    /**
     * @brief Initialize the lattice structure.
     * 
     * This virtual function must be implemented by derived classes
     * to set up the lattice according to specific rules or data.
     */

    virtual void initialize() = 0;

        /**
     * @brief Evaluate and return the energy of the lattice.
     * 
     * This function computes and returns the current energy state
     * of the lattice structure.
     * 
     * @return float The energy value of the lattice.
     */

    virtual float evaluate_energy() const = 0;

        /**
     * @brief Print the current state of the lattice.
     * 
     * This function provides a visualization or textual representation
     * of the lattice's current state.
     */

    virtual void print_lattice() const = 0;
        /**
     * @brief Get the interaction energy of the lattice elements.
     * 
     * This function calculates and returns the energy due to interactions
     * between elements of the lattice.
     * 
     * @return float The interaction energy of the lattice.
     */

    virtual float get_interaction_energy() const = 0; 

        /**
     * @brief Virtual destructor for the AbstractLattice class.
     * 
     * Ensures that derived class destructors are called correctly.
     */


    virtual ~AbstractLattice() = default;
    
   
};

#endif // ABSTRACTLATTICE_H