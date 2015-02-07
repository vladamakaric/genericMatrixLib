# genericMatrixLib
Experimental generic matrix/vector library in C++ using new C++11 features. The idea is to create a generic (arbitrary dimension and type) matrix/vector library in which dimensional mismatches would issue compile time errors. For example if you create a 7x2 matrix, the multiplication operator code is generated only for 2xN matrices, where N is any natural number. Also I began implementing implicit (upward) type conversion. For example if you multiply a 2x2 matrix of doubles with a 2x2 matrix of floats or ints, the result is a 2x2 double valued matrix, this is achieved with the new C++11/14 features. 
