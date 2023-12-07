#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include<fstream>
#include<ctime>
#include <chrono>
#include <omp.h>



#define L 100
#define N (L*L)
int A = 25; //Block size assigned to each thread
#define J 1.00
#define IT 6*1e7 //number of iterations

void print_lattice(std::vector < std::vector<int> >& matrix) {

    int i, j;
    for (i = 0; i < L; i++) {
        std::cout << std::endl;
        for (j = 0; j < L; j++) {
            if (matrix[i][j] == -1) {
                std::cout << "o" << " ";
            } else {
                std::cout << "x" << " ";

            }
        }
    }
}

//if needed function to correctly evaluate energy
float evaluate(std::vector<int>& lattice) {
    int sum=0;
    for (int i=0;i<N;i++){
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L]*lattice[i]*2;
        }
        if (i%L != 0) { //NO FIRST COLUMN
            sum += lattice[i-1]*lattice[i]*2;
        }

        if(i >= L*(L - 1)) { //LAST ROW
            sum += lattice[i-L*(L-1)]*lattice[i]*2; //times 2 two count also contribute where i = 0
        }
        if((i + 1) % L  == 0) { //LAST COLUMN
            sum += lattice[i-(L-1)]*lattice[i]*2;
        }
    }
    return -J*sum;
}

int flip(std::vector <int> & lattice, std::vector<float>& prob,float& energy, int n, bool color) {

    int sum = 0;
    int ID = omp_get_thread_num();
    int site = A * 2 * ID; //find from n the site to flip based on thread ID
    if(color){
        site +=( (ID/3) % 2) * A; //ID/3 gives the integer truncated part
    }
    else{
        site +=1 - ( (ID/3) % 2) * A;
    }

    if (site < L) {
        sum += lattice[site+L*(L-1)];
    }
    else {
        sum += lattice[site-L];
    }
    if (site % L == 0) {
        sum += lattice[site + (L - 1)];
    }
    else {
        sum += lattice[site - 1];
    }

    if (site >= L*(L - 1)) {
        sum += lattice[site - L*(L-1)];
    }
    else {
        sum += lattice[site + L];
    }
    if ((site+1) % L == 0) {
        sum += lattice[site - (L-1)];
    }
    else {
        sum += lattice[site + 1];
    }
    int delta = 2*sum*lattice[site];
    if (delta <= 0) {
        lattice[site] = -lattice[site];
    }
    else if (delta == 4) {
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[0] ){
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[1]) {
            lattice[n] = -lattice[n];
        }
        else{
            return 0;
        }
    }

    return 2*lattice[n];

}

void initialize_lattice(std::vector<int>& lattice, float& energy, int& M) {
    int k;
    int size = N;
    int sum = 0;
    for (int i = 0; i < size; i++) {
        k = rand() % 2;
        if (k == 0) {
            lattice[i] = -1;
            M -= 1;
        }
        else {
            lattice[i] = 1;
            M += 1;
        }
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L]*lattice[i]*2;
        }
        if (i%L != 0) { //NO FIRST COLUMN
            sum += lattice[i-1]*lattice[i]*2;
        }

        if(i >= L*(L - 1)) { //LAST ROW
            sum += lattice[i-L*(L-1)]*lattice[i]*2; //times 2 two count also contribute where i = 0
        }
        if((i + 1) % L  == 0) { //LAST COLUMN
            sum += lattice[i-(L-1)]*lattice[i]*2;
        }

    }
    energy+= -J*sum;
}






int write_file(std::vector<float>& energy_vec, std::vector<float>& m, std::vector<int>& t){
    std::ofstream myfile ("data.txt");
    if (myfile.is_open())
    {
        for(int i = 0; i < t.size(); i ++){
            myfile << energy_vec[i] << " "<<abs(m[i])<<' '<<t[i]<<'\n';
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";
    return 0;
}



float simulate(float T,std::vector <int> & lattice, float& energy, int& M,std::vector<int>& rand_vect, const int numThread) {

    using namespace std;
    vector<float> prob(2);
    prob[0] = exp(-4 * J/T);
    prob[1] = exp(-8 * J/T);

    vector<float> energy_vec(1);
    energy_vec[0] =  energy;

    vector<float> m(1);
    m[0] = (float)M/N;
    vector<int> t_axis(1);
    t_axis[0] = 0;
    #pragma omp parallel num_threads (numThread)
    {
        #pragma omp for reduction(+:M)
        for (int i = 0; i < IT-1; i+=2){
            int M_local = 0;
            M_local += flip(lattice, prob, energy, rand_vect[i], true); //true è nero
            //#pragma omp barrier
            M_local += flip(lattice, prob, energy, rand_vect[i+1], false); // false è bianco

            M += M_local;


        }
    }
    return (float)M/N;
}

void create_rand_vect(std::vector<int>& rand_vect_0) {
    int i;
    for (i = 0; i < IT; i++) {
        rand_vect_0.push_back(rand() % L );
    }
}

const int setBlockSize(int dimSideBlock){
    if(L%dimSideBlock == 0){
        const int numCellRow = L/dimSideBlock; // = numero di celle per colonna
        if(numCellRow % 2 == 0){
            return numCellRow*numCellRow/2;
        }
        else{
            std::cout<<"Con il valore inserito non posso definire una scacchiera";
        }
    }
    else{
        std::cout<<"La dimensione del reticolo non è multiplo della dimensione del blocco";
    }
    return 0;
}


int main() {
    using namespace std;
    unsigned seed = time(0);
    srand(seed);
    vector<int> lattice(N);
    const int numThread = setBlockSize(A);
    float energy = 0;
    int M =0;
    initialize_lattice(lattice,energy,M);
    float T = 0.1;
    vector<float> results (1);
    results[0] = 1;

    std::vector<int> rand_vect;
    create_rand_vect(rand_vect);

    while(T<=3.5){
        auto start = std::chrono::high_resolution_clock::now();  // Start timing before simulation
        results.push_back(simulate(T,lattice,energy,M,rand_vect,numThread));
        auto end = std::chrono::high_resolution_clock::now();  // End timing after simulation
        std::chrono::duration<double> elapsed = end - start;  // Calculate elapsed time

        T += 0.1;

        cout << "Time taken for simulation at temperature " << T << ": " << elapsed.count() << " seconds" << endl;
        cout << abs(results.back()) << endl;
    }



    return 0;

}

