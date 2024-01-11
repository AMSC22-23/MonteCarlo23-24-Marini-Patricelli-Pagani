// AbstractLattice.h

#ifndef ABSTRACTLATTICE_H
#define ABSTRACTLATTICE_H

class AbstractLattice {
public:
    virtual void initialize() = 0;
    virtual float evaluateEnergy() const = 0;
    virtual void printLattice() const = 0;
    virtual float getInteractionEnergy() const = 0; 
    virtual ~AbstractLattice() = default;
    
   
};

#endif // ABSTRACTLATTIC