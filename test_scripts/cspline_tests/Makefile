CXX				= g++
CXXFLAGS	= -g -O1 -Wall -Wextra
INC				= -I /usr/local/include -I /usr/local/lib/include \
						-I ../../src
LIBS			= -L/usr/local/lib -lm
PROGRAM		= test_cspline
OBJS			= test1.o cspline.o

$(PROGRAM): $(OBJS)
		$(CXX) $(CXXFLAGS) $(OBJS) -o $(PROGRAM) $(LIBS)

test1.o: test1.cpp ../../src/cspline.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c test1.cpp -o test1.o $(LIBS)

cspline.o: ../../src/cspline.cpp ../../src/cspline.hpp
	$(CXX) $(CXXFLAGS) $(INC) -c ../../src/cspline.cpp -o cspline.o $(LIBS)

.PHONY:
clean:
	/bin/rm -f $(OBJS)
	/bin/rm -f $(PROGRAM)
	/bin/rm -f core
