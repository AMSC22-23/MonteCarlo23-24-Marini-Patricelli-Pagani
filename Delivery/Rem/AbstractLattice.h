// AbstractLattice.h

#ifndef ABSTRACTLATTICE_H
#define ABSTRACTLATTICE_H

class AbstractLattice {
public:
    virtual void initialize() = 0;
    virtual float evaluate_energy() const = 0;
    virtual void print_lattice() const = 0;
    virtual float get_interaction_energy() const = 0; 
    virtual ~AbstractLattice() = default;
    
   
};

#endif // ABSTRACTLATTIC