CXX				= g++
CXXFLAGS	= -g -O1 -Wall -Wextra -std=c++11
INC				= -I /usr/local/include -I /usr/local/lib/include \
						-I ../../src
LIBS			= -L/usr/local/lib -lm
PROGRAM		= test_integrator
OBJS			= test1.o integrator.o

$(PROGRAM): $(OBJS)
		$(CXX) $(CXXFLAGS) $(OBJS) -o $(PROGRAM) $(LIBS)

test1.o: test1.cpp test1.hpp ../../src/integrator.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c test1.cpp -o test1.o $(LIBS)

integrator.o: ../../src/integrator.cpp ../../src/integrator.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c ../../src/integrator.cpp -o integrator.o $(LIBS)

.PHONY:
clean:
	/bin/rm -f $(OBJS)
	/bin/rm -f $(PROGRAM)
	/bin/rm -f core
