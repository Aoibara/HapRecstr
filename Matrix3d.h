#ifndef Matrix3d_h
#define Matrix3d_h

#include<map>
#include "Matrix615.h"
using namespace std;


template <class T>
class Matrix3d {
public:
    map< int, Matrix615<T> > mdata;
    Matrix3d(){
    } // flexibility to set a default matrix not specifying nrow ncol
    Matrix3d(int nrow, int ncol, int nslice, T val = 0) {
        Matrix615<T> temp(nrow, ncol, val); // make n rows
        for(int i=0; i < nslice; ++i) {
            mdata.insert(pair<int, Matrix615<T> > (i,temp)); // make n cols with default value val
        }
    }
    
    int sliNums() {return (int)mdata.size();}
};





#endif /* Matrix3d_h */
