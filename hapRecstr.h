#ifndef hapRecstr_h
#define hapRecstr_h

#include<set>
#include<algorithm>
#include<iomanip>
#include<utility>
#include "Matrix615.h"
#include "Matrix3d.h"
#include <boost/random/uniform_real.hpp> 
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#define TINY 1e-6
#define ZEPS 1e-10

using namespace std;

struct comparep {
    bool operator()(const pair<int,int> &left, const pair<int,int> &right) {
        return left.second > right.second;
    }
};

class HAPR {

public:
    // parameters
    int k ; // K : number of possible alleles per loci
    int n;  // n : number of observations
    int m ; // m : number of marker locis
    
    vector<double> pis; // 1*k. initial states
    Matrix615<int> g; // n*2m. observed genotypes
    Matrix615<double> trans; // (m-1)*(k*k). trans.data[0][i*k+j] corresponds to transmission from state[i]@loci[0] to state[j]@loci[1]
    Matrix615<double> emis; // m*(k*k). emis.data[0][i*k+j] corresponds to emission from state[i]@loci[0] of allele[j]
    
    // storages for dynamic programming
    Matrix615<double> alphas, betas, deltas; //  m*(k*k)
    Matrix615<double> probt, probe; // m*(k*k), (m-1)*(k*k). gammas. probability matrix for updating trans and emis, prob.data[0][i*k+j] corrensponds to the probability of getting state a=i and state/allele b=j @loci[0] 
    Matrix3d<int> phis;  // to backtrack optimal path
    Matrix615<int> paths;
    map< vector<int>, int> unihaps; // unique haplotypes and count

    // constructor: states, obs, markers
    HAPR(int states, int obs, int markers) : k(states), n(obs),
    m(markers), trans(markers-1, states*states , 0), emis(markers, states*states, 0), alphas(markers, states*states, 0), betas(markers, states*states, 0), deltas(markers, states*states, 0), probt(markers, states*states, 0), probe(markers, states*states, 0),phis(markers, states*states, 2 , -1), paths(markers, 2,-1){
        pis.resize(k); //
    }
    
    void initParams(double ro, double nu, double eta);
    
    void forward(int id);
    void backward(int id);
    void forwardBackward(int id); // E-step
    void updateParams(); //M-step including all above
    double locusLK(int id, int locus);
    double totalLK();
    int check_tol(double fmax, double fmin, double ftol); // check tolerance
    void runEM(double eps, int maxloop); // *public call function*
    void viterbi(int id);
    void hapOut(const char* fileName1, const char* fileName2);// *public call function*
    void sortOut2(vector< pair<int,int> > &pairs);
};

void HAPR::initParams(double ro, double nu, double eta){
    double t1 = 1 - ro, t2 = ro / (k-1);
    double e1 = 1 - nu, e2 = nu / (k);
    
    for(int i=0; i< k ; ++i){
        pis[i]= 1. / k ;
    }
    //initialize trans
    for(int t=0; t<m-1; ++t){
        for(int i=0; i<k; ++i ){
            for(int j=0; j<k; ++j){
                if(i==j){
                    trans.data[t][i*k+j] = t1;
                }
                else{
                    trans.data[t][i*k+j] =t2;
                }
            }
        }
    }
    
    //initialize emis, find and assign major alleles
    vector<int> permutation, best;
    double locallike =0, max =0;
    for(int i=0; i<k; ++i)
        permutation.push_back(i);
    best.resize(k);
    //t=0
    while(next_permutation(permutation.begin(), permutation.end())){
        for(int i=0; i<k; ++i ){
            for(int j=0; j<k; ++j){
                if(j==permutation[i]){
                    emis.data[0][i*k+j] = e1;
                }
                else{
                    emis.data[0][i*k+j] = e2;
                }
            }
        }
        for(int id=0; id<n; ++id){ // find local maximum
            int g1= g.data[id][0], g2= g.data[id][1];
            for(int i=0; i<k; ++i){ // state=i
                locallike += pis[i] * 0.5 * (emis.data[0][g1] + emis.data[0][g2]);
            }
        }
        if(locallike > max){
            max = locallike;
            best = permutation;
        }
    }
    for(int i=0; i<k; ++i ){ //assign major allele to states at current loci
        for(int j=0; j<k; ++j){
            if(j==best[i]){
                emis.data[0][i*k+j] = e1;
            }
            else{
                emis.data[0][i*k+j] = e2;
            }
        }
    }
    //t>=1
    for(int t=1; t<m; ++t){
        while(next_permutation(permutation.begin(), permutation.end())){
            locallike =0;
            for(int i=0; i<k; ++i ){
                for(int j=0; j<k; ++j){
                    if(j==permutation[i]){
                        emis.data[t][i*k+j] = e1;
                    }
                    else{
                        emis.data[t][i*k+j] = e2;
                    }
                }
            }
            
            // get local max likelihood
            for(int id=0; id<n; ++id){ // find local maximum
                int g1= g.data[id][2*t], g2= g.data[id][2*t+1];
                for(int i=0; i<k; ++i){ // state=i
                    int temptrans = 0;
                    for(int j=0; j<k; ++j){ // columnwise sum of trans prob to current state(i)
                        temptrans += trans.data[t-1][j*k+i];
                    }
                    locallike += temptrans * 0.5 * (emis.data[t][g1] + emis.data[t][g2]);
                }
            }
            if(locallike > max){
                max = locallike;
                best = permutation;
            }
        }
        for(int i=0; i<k; ++i ){ //assign major allele to states at current loci
            for(int j=0; j<k; ++j){
                if(j==best[i]){
                    emis.data[t][i*k+j] = e1;
                }
                else{
                    emis.data[t][i*k+j] = e2;
                }
            }
        }
    }
    
    // perturb
    mt19937 rng;
    rng.seed(time(0));
    uniform_real<> uni(-eta,eta);
    double rowsum =0;
    for(int t=0; t<m; ++t){
        for(int i=0; i<k; ++i){
            for(int j=0; j<k; ++j){
                emis.data[t][i*k+j] *= exp(uni(rng)); // *exp(X)
                rowsum += emis.data[t][i*k+j];
            }
            for(int j=0; j<k; ++j){
                emis.data[t][i*k+j] /= rowsum; // normalize
            }
            rowsum=0;
        }
    }
    for(int t=0; t<m-1; ++t){
        for(int i=0; i<k; ++i){
            for(int j=0; j<k; ++j){
                trans.data[t][i*k+j] *= exp(uni(rng)); // *exp(X)
                rowsum += trans.data[t][i*k+j];
            }
            for(int j=0; j<k; ++j){
                trans.data[t][i*k+j] /= rowsum; // normalize
            }
            rowsum =0;
        }
    }
    rowsum=0;
    for(int i=0; i<k; ++i){
        pis[i] *= exp(uni(rng));
        rowsum += pis[i];
    }
    for(int i=0; i<k; ++i){
        pis[i] /= rowsum;
    }
}


void HAPR::forward(int id=0) { //left probability of emitting the initial segment g~t-1 and ending at at,bt. alphas=L
    Matrix615<double> sum1(m,k*k,0);// see article 3.3
    int g1= g.data[id][0], g2= g.data[id][1];
    for(int i=0; i < k; ++i) { // t=0, a(t)=i
        for(int j=0; j<k; ++j){ // t=0, b(t)=j
            alphas.data[0][i*k+j] = pis[i] * pis[j];
            for(int l=0; l<k; ++l){ // a(t)=i, b(t)=j, b(t+1)=l
                sum1.data[0][i*k+l] += alphas.data[0][i*k+j] * trans.data[0][j*k+l] * 0.5 * (emis.data[0][i*k+g1] * emis.data[0][j*k+ g2] + emis.data[0][i*k+ g2] * emis.data[0][j*k+ g1]);
            }
            
        }
    }
    for(int t=1; t < m; ++t) { //t>0
        g1= g.data[id][2*t], g2= g.data[id][2*t+1];
        for(int i=0; i < k; ++i) { //a(t)=i
            for(int j=0; j < k; ++j) { //b(t)=j
                alphas.data[t][i*k+j] = 0;
                for(int l=0; l<k; ++l){ //a(t-1)=l, b(t+1)=l
                    alphas.data[t][i*k+j] += trans.data[t-1][l*k+i] * sum1.data[t-1][l*k+j];
                    if(t<m-1){
                        sum1.data[t][i*k+l] += alphas.data[t][i*k+j] * trans.data[t][j*k+l] * 0.5 * (emis.data[t][i*k+g1] * emis.data[t][j*k+g2] + emis.data[t][i*k+g2] * emis.data[t][j*k+g1]);
                    }
                }
            }
        }
    }
}


void HAPR::backward(int id=0) { //right probability of emitting end segment g(t+1)~g(m) from state a(t),b(t). betas=R
    Matrix615<double> sum2(m,k*k,0);
    int g1=g.data[id][2*m-2], g2=g.data[id][2*m-1];
    for(int i=0; i < k; ++i) { //t=m-1, a(t)=i
        for(int j=0; j<k; ++j){ // t=m-1, b(t)=j
            betas.data[m-1][i*k+j] = 1;
            for(int l=0; l<k; ++l){ // a(t)=i, b(t)=j, b(t-1)=l
                sum2.data[m-1][i*k+l] += betas.data[m-1][i*k+j] * trans.data[m-2][l*k+j] * 0.5 * (emis.data[m-1][i*k+g1] * emis.data[m-1][j*k+ g2] + emis.data[m-1][i*k+g2] * emis.data[m-1][j*k+g1]);
            }

        }
        
    }
    for(int t=m-2; t >=0; --t) {
        int g1= g.data[id][2*t], g2= g.data[id][2*t+1];
        for(int i=0; i < k; ++i) { //a(t)=i
            for(int j=0; j < k; ++j) { //b(t)=j
                betas.data[t][i*k+j] = 0;
                for(int l=0; l<k; ++l){ //a(t+1)=l, b(t-1)=l
                    betas.data[t][i*k+j] += trans.data[t][i*k+l] * sum2.data[t+1][l*k+j];
                    if(t>0){
                    sum2.data[t][i*k+l] += betas.data[t][i*k+j] * trans.data[t-1][l*k+j] * 0.5 * (emis.data[t][i*k+g1] * emis.data[t][j*k+g2] + emis.data[t][i*k+g2] * emis.data[t][j*k+g1]);
                    }
                }
            }
        }
    }
}


void HAPR::forwardBackward(int id=0) { //notmalization procedure
    forward(id);
    backward(id);
    
    
    // go through locus and calculate probe probt
    for(int t=0; t < m; ++t) {
        int g1= g.data[id][2*t], g2= g.data[id][2*t+1];
        double marginal = locusLK(id, t);
        
        // probe: unnormalized prob of state b emitting allele y @ locus t (numerator of evaluation formula 3, update value emis)
        for(int i=0; i<k; ++i){ //b=i
            for(int j=0; j<k; ++j){ //y=j
                double jointe = 0;
                for(int l=0; l<k; ++l){ // sum over a=l
                    jointe += alphas.data[t][l*k+i] * betas.data[t][l*k+i] * 0.5 * (emis.data[t][l*k+g1] * emis.data[t][i*k+g2] + emis.data[t][l*k+g2] * emis.data[t][i*k+g1]);
                }
                
                //cout << "jointe: " << jointe << ", "  << "marginal: " << marginal << ", "<< "j=" << j << ", "<< "i=" << i << ", "<< "t=" << t << endl;
                probe.data[t][i*k+j] = jointe / marginal;
            }
        }
        
        // probt:  unnormalized prob of state a@locus t transitting to state b@locus t+1 (evaluation formula 2, update value trans)
        if(t<m-1){
            for(int i=0; i<k; ++i){ //a=i @t
                for(int j=0; j<k; ++j){ //b=j @t+1
                    double jointt = 0;
                    for(int l=0; l<k; ++l){ // b=l @t
                        jointt += alphas.data[t][i*k+l] * betas.data[t][i*k+l] * trans.data[t][i*k+j];
                    }
                    probt.data[t][i*k+j] = jointt / marginal;
                }
            }
        }
    }
}


void HAPR::updateParams(){ //run through all data get new emis and trans for updating parameters
    Matrix615<double> newemis(m,k*k,0), newtrans(m,k*k,0);
    vector<double> sume(m), sumt(m); //keep track of sum of all prob for normalizing
    for(int id =0; id< n; ++id){
        forwardBackward(id); // update alphas, betas, probe, probt for g[id]
        for(int t=0; t<m; ++t){
            for(int i=0; i<k; ++i){
                for(int j=0; j<k; ++j){
                    newemis.data[t][i*k+j] += probe.data[t][i*k+j];
                    sume[t] += probe.data[t][i*k+j];
                    if(t<m-1){
                        newtrans.data[t][i*k+j] += probt.data[t][i*k+j];
                        sumt[t] += probt.data[t][i*k+j];
                    }
                }
            }
        }
    }
    
    // normalize newemis, newtrans
    // update emis, trans with newemis, newtrans
    for(int t=0; t<m; ++t){
        for(int i=0; i<k; ++i){
            for(int j=0; j<k; ++j){
                emis.data[t][i*k+j] = newemis.data[t][i*k+j] / sume[t]; //update emis
                if(t<m-1){
                    trans.data[t][i*k+j] = newtrans.data[t][i*k+j] / sumt[t]; // update trans
                }
            }
        }
    }
}

double HAPR::locusLK(int id, int locus) {
    double llk = 0.0;
    int g1= g.data[id][2*locus], g2= g.data[id][2*locus+1];
    for(int i=0; i<k; ++i){ //t=0, a=i, b=j
        for(int j=0; j<k; ++j){
            llk += alphas.data[locus][i*k+j] * 0.5* (emis.data[locus][i*k+g1] * emis.data[locus][j*k+ g2] + emis.data[locus][i*k+g2] * emis.data[locus][j*k+g1]);
        }
    }
    return llk;
}

double HAPR::totalLK() {
    double lk =1.0;
    for(int id=0; id<n; ++id){
        double glk = 0;
        for(int t=0; t<m; ++t){
            glk += locusLK(id,t);
        }
        lk *= glk;
    }
    return lk;
}

int HAPR::check_tol(double fmax, double fmin, double ftol) {
    double delta = fabs(fmax - fmin);
    double accuracy = (fabs(fmax) + fabs(fmin)) * ftol;
    return (delta < (accuracy + ZEPS));
}


void HAPR::runEM(double eps, int maxloop) {
    double lk = 0, prevLK = 0;
    int loop = 0;
    while( ( lk == 0 ) || ( check_tol(lk, prevLK, eps) == 0 ) ) {
        if(loop == maxloop) break;
        updateParams();
        prevLK = lk;
        lk = totalLK();
        loop++;
    }
}




void HAPR::viterbi(int id=0) {
    int g1= g.data[id][0], g2= g.data[id][1];
    // t=0, initialize delta matrix@loci1, delta[0*k+i]: emit g1 from state i, delta[1*k+i] emit g2
    for(int i=0; i < k; ++i) { // a=i
        for(int j=0; j<k; ++j){ //b=j
            deltas.data[0][i*k+j] = pis[i] * pis[j] * 0.5 * (emis.data[0][i*k+g1] * emis.data[0][j*k+ g2] + emis.data[0][i*k+ g2] * emis.data[0][j*k+ g1]);
        }
    }
    
    for(int t=1; t < m; ++t) { //loci=t>0
        g1= g.data[id][2*t], g2= g.data[id][2*t+1];
        for(int i=0; i < k; ++i) { // state a=i
            for(int j=1; j < k; ++j) { //state b=j
                int maxIda=0, maxIdb=0;
                double maxVal = 0;
                for(int l1=0; l1<k; ++l1){ //a'=l1
                    for(int l2=0; l2<k; ++l2){ //b'=l2
                        double val= deltas.data[t-1][l1*k+l2] * trans.data[t-1][l1*k+i] * trans.data[t-1][l2*k+j] * 0.5 * (emis.data[t][i*k+g1] * emis.data[t][j*k+ g2] + emis.data[t][i*k+ g2] * emis.data[t][j*k+ g1]);
                        if ( val > maxVal ) {
                            maxIda = l1;
                            maxIdb = l2;
                            maxVal = val;
                        }
                    }
                }
                deltas.data[t][i*k+j] = maxVal; //save delta
                phis.mdata[0].data[t][i*k+j] = maxIda;
                phis.mdata[1].data[t][i*k+j] = maxIdb; //save phi
            }
        }
    }
    
    
    // backtrack viterbi path
    paths.data[m-1][0] = 0;
    paths.data[m-1][1] = 0;
    double maxDelta = 0;
    for(int i=0; i < k; ++i) { //t=m-1
        for(int j=0; j < k; ++j){
            double deltaVal = deltas.data[m-1][i*k+j];
            if ( maxDelta < deltaVal ) {
                maxDelta = deltaVal;
                paths.data[m-1][0] = i;
                paths.data[m-1][1] = j;
            }
        }
    }
    for(int t=m-2; t >= 0; --t) { //t<m-1
        int ida = paths.data[t+1][0], idb = paths.data[t+1][1];
        paths.data[t][0] = phis.mdata[0].data[t+1][ ida*k+idb ];
        paths.data[t][1] = phis.mdata[1].data[t+1][ ida*k+idb ];
    }
}



/***dumped*** reconstruct haplotype for one specified id
void HAPR::hapGen(int id=0) {
    viterbi(id);
    for(int t=0; t<m; t++){
        int g1= g.data[id][2*t], g2= g.data[id][2*t+1];
        if(g1 == g2){ //if loci is homozygous
            hap.data[0][t]=g1;
            hap.data[1][t]=g2;
        }
        else{//if loci is heterozygous
            double p1emisg1 =emis.data[t][ paths.data[t][0] + g1];
            double p1emisg2 =emis.data[t][ paths.data[t][0] + g2];
            double p2emisg1 =emis.data[t][ paths.data[t][1] + g1];
            double p2emisg2 =emis.data[t][ paths.data[t][1] + g2];
            if(  p1emisg1*p2emisg2 > p1emisg2*p2emisg1  ){
                hap.data[0][t]=g1;
                hap.data[1][t]=g2;
            }
            else{
                hap.data[0][t]=g2;
                hap.data[1][t]=g1;
            }
        }
    }
    cout << "Original genotype for id " << id << ":" << endl;
    for(int i=0; i<2*m; ++i)
        cout << g.data[id][i] << " ";
    cout << endl;
    
    cout << "Reconstructed haplotype(s) for id " << id << ":" << endl;
    hap.print();
}
*/


// store reconstructed data into output file2
// fileName1: file to store reconstructed haplotype data, first column is id number, followed by haplotype (each id generate two haplotypes)
// fileName2: file to store summary of haplotype data, last column is count of the haplotype in the data population
void HAPR::hapOut(const char* fileName1, const char* fileName2){
    Matrix615<int> hap(2,m);
    vector< pair<int,int> > unihapp;
    ofstream myfile1 (fileName1);
    if (myfile1.is_open()){
        for(int id=0; id<n; ++id){
            viterbi(id);
            for(int t=0; t<m; t++){
                int g1= g.data[id][2*t], g2= g.data[id][2*t+1];
                if(g1 == g2){ //if loci is homozygous
                    hap.data[0][t]=g1;
                    hap.data[1][t]=g2;
                }
                else{//if loci is heterozygous
                    double p1emisg1 =emis.data[t][ paths.data[t][0] + g1];
                    double p1emisg2 =emis.data[t][ paths.data[t][0] + g2];
                    double p2emisg1 =emis.data[t][ paths.data[t][1] + g1];
                    double p2emisg2 =emis.data[t][ paths.data[t][1] + g2];
                    if(  p1emisg1*p2emisg2 > p1emisg2*p2emisg1  ){
                        hap.data[0][t]=g1;
                        hap.data[1][t]=g2;
                    }
                    else{
                        hap.data[0][t]=g2;
                        hap.data[1][t]=g1;
                    }
                }
            }
            for(int i=0; i < 2;i++){
                myfile1 << id+1 << " ";
                for(int j=0;j < m; j++){
                    myfile1 << hap.data[i][j] << " ";
                }
                ++unihaps[hap.data[i]];
                myfile1 << "\n";
            }
        }
        myfile1.close();
        cout << "Output 1 saved to " << fileName1 << endl;
    }
    else cout << "Unable to open file" << fileName1 << endl;
    
    sortOut2(unihapp);
    ofstream myfile2 (fileName2);
    if (myfile2.is_open()){
        map< vector<int>, int >::const_iterator it = unihaps.begin();
        for(int i=0; i<(int)unihapp.size(); ++i){
            advance(it,unihapp[i].first);
            vector<int> temp= it -> first;
            for(int j=0; j< (int)temp.size(); ++j){
                myfile2 << temp[j] << " ";
            }
            advance(it, 0- unihapp[i].first);
            myfile2 << unihapp[i].second / (2.*n) << "\n";
        }
        myfile2.close();
        cout << "Output 2 saved to " << fileName2 << endl;
    }
    else cout << "Unable to open file" << fileName2 << endl;
    
    // print parameters
    ofstream myfile3 ("emis.txt");
    if (myfile3.is_open()){
        for(int i=0; i<m; ++i){
            for(int j=0; j<k*k; ++j){
                myfile3 << emis.data[i][j] << " ";
            }
            myfile3 << endl;
        }
    }
    myfile3.close();
    ofstream myfile4 ("trans.txt");
    if (myfile4.is_open()){
        for(int i=0; i<m-1; ++i){
            for(int j=0; j<k*k; ++j){
                myfile4 << trans.data[i][j] << " ";
            }
            myfile4 << endl;
        }
    }
    myfile4.close();
    ofstream myfile5 ("pis.txt");
    if (myfile5.is_open()){
        for(int i=0; i<k; ++i){
            myfile5 << pis[i] << " ";
        }
    }
    myfile5.close();
}

void HAPR::sortOut2(vector< pair<int,int> > &pairs){
    pair<int,int> tokp;  //pair<unihapid, count>
    int toki=0;
    for(map< vector<int>, int >::const_iterator it = unihaps.begin(); it != unihaps.end(); ++it){
        tokp.first = toki;
        tokp.second = it -> second;
        pairs.push_back(tokp);
        ++toki;
    }
    sort(pairs.begin(),pairs.end(),comparep());
}




#endif /* hapRecstr_h */
