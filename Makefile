LIB=-L/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/  -L/opt/ibm/ILOG/CPLEX_Studio2211/concert/lib/x86-64_linux/static_pic/
INC=-I/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio2211/concert/include

CFLAGS=-Wall $(LIB) $(INC) -std=c++0x -m64 -g  -fno-strict-aliasing #-O3
CPPFLAGS=-DNDEBUG -DIL_STD
CXX=g++
main: main.o BFieldThick.o   BFieldThick_Matrix.o Integral.o Ellipse.o Model2.o Model.o supplementary_functions.o
	$(CXX) $(CFLAGS) -o main main.o BFieldThick.o  BFieldThick_Matrix.o Integral.o Ellipse.o Model2.o Model.o supplementary_functions.o -lilocplex -lconcert -lcplex -lm -lpthread

main.o: main.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c main.cpp
BFieldThick.o: BFieldThick.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c BFieldThick.cpp
BFieldThick_Matrix.o:  BFieldThick_Matrix.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c  BFieldThick_Matrix.cpp
Integral.o: Integral.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c Integral.cpp
Ellipse.o: Ellipse.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c Ellipse.cpp
Model2.o: Model2.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c Model2.cpp
Model.o: Model.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c Model.cpp
supplementary_functions.o: supplementary_functions.cpp
	g++ $(CFLAGS) $(CPPFLAGS) -c supplementary_functions.cpp

clean:
	rm -f main main.o BFieldThick.o  BFieldThick_Matrix.o  Integral.o Ellipse.o Model2.o Model.o supplementary_functions.o
