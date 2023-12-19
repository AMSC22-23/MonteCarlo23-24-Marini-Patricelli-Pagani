# General notes
- I am focusing just of the **Hands-on** folder
- In the readme I would like to have few examples of how to run the code and what should happen
- Makefile must have the correct capitalization in Linux and does not work

# Code
## Major
- You still have not removed the `muparser-2.3.4` folder even if now there is the submodule
- You have not deleted the `HyperRectangle.cpp` and `HyperShpere.cpp` files even if they are empty
- Since `HyperRectangle.hpp` and `HyperShpere.hpp` are not template, the definition of methods should be in a source
- There is a missing `#` at line 4 of `HyperRectangle.hpp`
- You should remove the `main` executable from the folder
- I suggest to avoid mixing English and Italian in naming methods and attributes
- Things that could be passed by reference are passed by copy (in partigulare `std::string`)
- Number of OMP threads should not be hard-coded
- At line 57 of `HyperShpere.hpp`, instead of creating at each iteration the engine with the custom range, you could create it just once with distribution in [0, 1] and then rescale the random number in the correct interval by means of a linear transformation
- At line 54 of `HyperShpere.hpp`, the parallelization is counterproductive since the function `generateVector` is already called in a parallel section, so you are splitting again resources even if they are not physically avaiable
## Minor
- The `.at()` method of arrays is a bit slower than the `operator[]` since has bound checks, prefer the latter for vectors. If you want to have memory safery use the `-fsanitize=address` at compile time
- Parameters file should not be in the source code folder, but in a separate folder