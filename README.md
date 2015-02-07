# genericMatrixLib
Experimental generic matrix/vector library in C++ using new C++11 features. The idea is to create a generic matrix/vector library in which dimensional mismatches would issue compile time errors. For example if you create a 7x2 matrix, the multiplication operator code is generated only for 2xN matrices, where N is any natural nuber. 
