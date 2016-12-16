#ifndef _Matrix_615_H_ //avoid double compile
#define _Matrix_615_H_

#include<iostream>
#include<vector>
#include<string>
#include<cstdlib>
#include<climits>
#include<cmath>
#include<boost/tokenizer.hpp>
#include<fstream>
#include<boost/lexical_cast.hpp>
using namespace std;
using namespace boost;

template <class T>
class Matrix615 {
public:
    vector<vector<T> > data;
    Matrix615(){
    } // flexibility to set a default matrix not specifying nrow ncol
    Matrix615(int nrow, int ncol, T val = 0) {
        data.resize(nrow); // make n rows
        for(int i=0; i < nrow; ++i) {
            data[i].resize(ncol,val); // make n cols with default value val
        }
    }
    
    int rowNums() {return (int)data.size();}
    int colNums() {return (data.size()==0) ? 0: (int)data[0].size();}
    void readFromFile(const char* fileName);
    void readFromFileint(const char* fileName);
    void print();
};


template <class T>
void Matrix615<T>::print(){
    for(int i=0; i < rowNums();i++){
        for(int j=0;j < colNums(); j++){
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

template <class T>
void Matrix615<T>::readFromFile(const char* fileName) {
    ifstream ifs(fileName);
    if ( ! ifs.is_open() ) {
        cerr<< "Cannot open file" <<fileName<< endl;
        abort();
    }
    string line;
    char_separator<char> sep(", \t \"");
    typedef tokenizer<char_separator<char> > wsTokenizer;
    data.clear();
    int nr=0, nc=0;
    while(getline(ifs, line) ) {
        if (line[0]=='#') continue;
        wsTokenizer t(line,sep);
        data.resize(nr+1);
        for(wsTokenizer::iterator i=t.begin(); i !=t.end(); ++i) {
            data[nr].push_back(lexical_cast<T>(i->c_str()));
            if (nr==0) ++nc;
        }
        if (nc != (int)data[nr].size() ) {
            cerr<<"The input file is not rectangle at line "<<nr<<endl;
            abort();
        }
        ++nr;
    }
}



// further check: integer, column number even
template <class T>
void Matrix615<T>::readFromFileint(const char* fileName) {
    
    ifstream ifs(fileName);
    if ( ! ifs.is_open() ) {
        cerr<< "Cannot open file" <<fileName<< endl;
        abort();
    }
    string line;
    char_separator<char> sep(", \t \"");
    typedef tokenizer<char_separator<char> > wsTokenizer;
    data.clear();
    int nr=0, nc=0;
    while(getline(ifs, line) ) {
        if (line[0]=='#') continue;
        wsTokenizer t(line,sep);
        data.resize(nr+1);
        for(wsTokenizer::iterator i=t.begin(); i !=t.end(); ++i) {
            try{
                data[nr].push_back(lexical_cast<int>(*i));
            }
            catch(const bad_lexical_cast &){
                cerr<<"Elements in the input file should be integers: line " << nr << " column "<< nc <<endl;
                abort();
            };
            if (nr==0) ++nc;
        }
        if (nc != (int)data[nr].size() ) {
            cerr<<"The input file is not rectangle at line "<<nr<<endl;
            abort();
        }
        if (nc % 2 !=0){
            cerr<<"Input file should have even number of columns "<<nr<<endl;
            abort();
        }
        ++nr;
    }
}

#endif
