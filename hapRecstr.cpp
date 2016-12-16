#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include "Matrix615.h"
#include "Matrix3d.h"
#include "hapRecstr.h"

using namespace std;

int main(int argc, char** argv){
    /*
    Matrix3d<int> test2(2,2,3);
    for(int i=0; i<2;i++){
        test2.mdata[i].print();
        cout << endl;
        cout << endl;
    }
     */
    if(argc<=4){
        cerr<< "The program takes at least 4 arguments\n arg1: source file name (genotype data)\n arg2: output file 1 name (reconstructed haplotype data)\n arg3: output file 2 name (summary of unique haplotypes and frequency)\n arg4: number of states per SNP (entries in genotype data should not exceed this number)\n arg5 (optional): file containing emission\n arg6 (optional): file containing transmision\n arg7 (optional): file containing pis (first transmission)\n" << endl;
        abort();
    }
    
    Matrix615<int> input;
    // read from file and get k,n,m
    input.readFromFileint(argv[1]);
    int n=input.rowNums(), m=input.colNums()/2, k=(atoi)(argv[4]);
    if(argc>5 && argc<8){
        cerr << "3 files for parameter initialization required" << endl;
        abort();
    }
    if(argc>8){
        cerr << "The program takes at most 7 arguments" << endl;
        abort();
    }
    
    // pass to HAPR object
    HAPR sample(k,n,m);
    sample.g= input;
    
    // parameter from file
    if(argc==8){
        sample.emis.readFromFile(argv[5]);
        sample.trans.readFromFile(argv[6]);
        
        Matrix615<double> temp;
        temp.readFromFile(argv[7]);
        for(int i=0; i< temp.colNums(); ++i){
            sample.pis[i] = temp.data[0][i];
        }
    }
    
    // parameter default, user specify constants ro, nu, eta
    else{
        double ro=-1., nu=-1.,eta=-1.;
        cout << "Set ro, nu, eta for parameter initialization (type 0 to use default value ro=0.1, nu=0.01, eta=0.8):" << endl;
        cin >> ro >> nu >> eta;
        while(cin.fail()){
            cin.clear();
            cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            cerr << "Input not valid\n Set ro, nu, eta for parameter initialization (type 0 to use default value ro=0.1, nu=0.01, eta=0.8):";
            cin >> ro >> nu >> eta;
        }
        while(ro < 0 || ro >=1 || nu < 0 || nu >=1 || eta < 0 || eta >=1){
            cerr << "Input not valid\n Set ro, nu, eta for parameter initialization (type 0 to use default value ro=0.1, nu=0.01, eta=0.8):" << endl;
            cin >> ro >> nu >> eta;
        }
        if(ro==0)
            ro=0.1;
        if(nu==0)
            nu= 0.01;
        if(eta==0)
            eta=0.8;
        sample.initParams(ro,nu,eta);
    }
    
    // set convergence threshold
    double in1= -1.;
    cout << "Model initialized\n Set convergence threshold (type 0 to use default value 1e-6):" << endl;
    cin >> in1;
    while(cin.fail()){
        cin.clear();
        cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        cerr << "Input not valid\n Set convergence threshold (type 0 to use default value 1e-6):";
        cin >> in1;
    }
    while(in1 < 0 || in1 >=1){
        cerr << "Input not valid\n Set convergence threshold (type 0 to use default value 1e-6):" << endl;
        cin >> in1;
    }
    if(in1==0)
        in1=1e-6;
    
    // set maximum EM loops
    int in2= -1;
    cout << " Set maximum EM restarts (type 0 to use default value 25):" << endl;
    cin >> in2;
    while(cin.fail()){
        cin.clear();
        cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        cerr << "Input not valid\n Set maximum EM loops (type 0 to use default value 25):";
        cin >> in2;
    }
    while(in2 < 0){
        cerr << "Input not valid\n Set maximum EM loops (type 0 to use default value 25):" << endl;
        cin >> in2;
    }
    if(in2==0)
        in2=25;
    
    
    sample.runEM(in1, in2);
    sample.hapOut(argv[2], argv[3]);


    
    
    return 0;
}
