CXX				= g++
CXXFLAGS	= -g -o3 -Wall -Wextra -std=c++11 -fopenmp
INC				= -I/usr/local/include -I/usr/local/lib/include
LIBS			= -L/usr/local/lib -lm
PROGRAM		= run
OBJS			= utils.o integrator.o therm_quantities.o potentials.o\
						config_integrals.o cspline.o eam.o simulate.o

$(PROGRAM): $(OBJS)
		$(CXX) $(CXXFLAGS) $(OBJS) -o $(PROGRAM) $(LIBS)

utils.o: utils.cpp utils.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c utils.cpp -o utils.o $(LIBS)

integrator.o: integrator.cpp integrator.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c integrator.cpp -o integrator.o $(LIBS)

therm_quantities.o: therm_quantities.cpp therm_quantities.hpp utils.hpp\
	                  config_integrals.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c therm_quantities.cpp -o therm_quantities.o \
  $(LIBS)

potentials.o: potentials.cpp potentials.hpp utils.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c  potentials.cpp -o potentials.o $(LIBS)

config_integrals.o: config_integrals.cpp config_integrals.hpp utils.hpp\
	                  integrator.hpp potentials.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c config_integrals.cpp -o config_integrals.o \
		                 $(LIBS)

cspline.o: cspline.cpp cspline.hpp 
	$(CXX) $(CXXFLAGS) $(INC) -c cspline.cpp -o cspline.o $(LIBS)

eam.o: eam.cpp eam.hpp utils.hpp cspline.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c eam.cpp -o eam.o $(LIBS)

simulate.o: simulate.cpp simulate.hpp utils.hpp therm_quantities.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c simulate.cpp -o simulate.o $(LIBS)

.PHONY:
clean:
	/bin/rm -f $(OBJS)
	/bin/rm -f $(PROGRAM)
	/bin/rm -f core
