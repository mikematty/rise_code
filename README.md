This repository contains code for computing various thermodynamic
quantities and other information relating to the study of close packing solids
and vacancies in close packing solids.  The main utility is the Stillinger
expansion of the canonical partition function.

REQUIREMENTS:
  Open MP
  C++ compiler that is capable of compiling using the C++ 11 standard
    (preferably, g++ 4.7.3 as this is what this software is tested on)

GETTING THE RIGHT COMPILER:
  If you don't have the right version of gcc/g++, getting it and making it the
    default can be a real pain.  Here is a line-by-line of what you need to do:

  sudo add-apt-repository ppa:ubuntu-toolchain-r/test
  sudo apt-get update
  sudo apt-get install gcc-4.7 g++-4.7
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.7 50
  sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.7 50

  To make sure it worked, typing g++ --version should print 4.7.3 at the end
    of the first line of output.

PARTICLE INDEXING:
  One of the most important parts of this code is the indexing scheme for the
    particles in the (fcc) lattice.  The drawing below illustrates the
    in-plane nearest neighbors of a particle.  all particles shown below are
    indexed by (x,y), and have z-coordinate 0
                        ___
                    ___/10 \___
            -1,-1->/   \___/-11\   x
                   \___/00 \___/   |_y
                   /0-1\___/01 \    
                   \___/-10\___/
                       \___/
    
  Out of plane neighbors are such that (0,0,1) is centered on the vertex
    of particles (0,0,0),(1,0,0), and (-1,1,0), while (0,0,-1) is centered
    on particles (0,0,0),(-1,0,0), and (1,-1,0).  The in-plane indexing
    is the same across all planes, while the position of the (0,0,n) particle
    for plane with z-index n is determined by the congruence class of n
    modulo 3, with the different cases as illustrated by the 0, +/-1 cases
    given above.

POTENTIAL CLASSES:
  A complete discussion of the implementation of various potentials is
    slightly premature at this point, but this is an appropriate point to
    discuss the extension of the software to include new potentials.
    In the file potentials.hpp, one can find a declaration of an abstract
    base class called Potential.  This class has two pure virtual functions
    with the following call signatures:
      double exptential_1p(vector<int*> cs, sim_i *sim, double *coords);
      double exptential_2p(vector<int*> cs1, vector<int*> cs2,
                           sim_i *sim, int *from2, double *coords);
    The detailed specs can be found in potentials.hpp.  The declaration of
    Potential also has a field called max_dist, which is a value of type
    double corresponding to the cutoff range of the potential.
    
  Any class derived from Potential with implemented exptential functions
    that sets a value for max_dist can be used as a potential for any of the
    simulations this software is capable of running.  Users should feel free
    to implement any additional potential classes as needed, although it is
    recommended to create separate files for them and merely #include the
    potentials.hpp header.

USER VARIABLES:
            
DISCUSSION OF FILES:
  This section contains a brief description of the utility of each function
    in each file.  For more detailed information such as calling signatures,
    argument descriptions, and return value discussions, consult the header
    file associated with a given function.

COMPILING:
  From the root directory of the repository, cd into 'src'.
  From src, run the make command
  To delete resultant object files and executable, run make clean

RUNNING:
  Compilation produces an executable called run
  All available simulations are executed by performing ./run T rho
  where T is an input temperature value, and rho is an input density value

CHANGING THE SIMULATION:
  Uh. we'll work this out
