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
#define IT 6*1e7//number of iterations

void print_lattice(std::vector <int> & lattice) {

    int i;
    for (i = 0; i < N; i++) {
            if(i%L == 0) std::cout<<std::endl;
            if (lattice[i] == -1) {
                std::cout << "o" << " ";
            } else {
                std::cout << "x" << " ";


            }
        }
    std::cout<<std::endl;
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
    int site = A * 2 * ID + n; //find from n and Thread ID the site to flip
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
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    return 2*lattice[site];

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
    prob[0] = exp(-4 * J / T);
    prob[1] = exp(-8 * J / T);


    vector<int> t_axis(1);
    t_axis[0] = 0;
#pragma omp parallel num_threads (numThread) reduction(+ : M)
    {
#pragma omp single
        {
            for (int i = 0; i < IT / (2 * numThread) - 1; i++) {
#pragma omp taskgroup task_reduction(+: M)
                {
                    for (int j = 0; j < 2 * numThread; j += 2) {
                        int index = i * 2 * numThread + j;
#pragma omp task in_reduction(+: M)
                        {
                            M += flip(lattice, prob, energy, rand_vect[index], true); //true è nero
                            M += flip(lattice, prob, energy, rand_vect[index + 1], false); //false è bianco
                        }
                    }
                }
            }
        }
    }
/*
     #pragma omp parallel for num_threads (numThread) reduction(+ : M)
     {
         for (int i = 0; i < IT-1; i+=2){
             int M_local = 0;
             M += flip(lattice, prob, energy, rand_vect[i], true); //true è nero
             M += flip(lattice, prob, energy, rand_vect[i+1], false); //true è nero

         }
     }*/


    float m = (float)M/N;
    return m ;
}

    void create_rand_vect(std::vector<int> &rand_vect_0) {
        int i;
        for (i = 0; i < IT; i++) {
            rand_vect_0.push_back(rand() % L);
        }
    }

    const int setBlockSize(int dimSideBlock) {
        if (L % dimSideBlock == 0) {
            const int numCellRow = L / dimSideBlock; // = numero di celle per colonna
            if (numCellRow % 2 == 0) {
                return numCellRow * numCellRow / 2;
            } else {
                std::cout << "Con il valore inserito non posso definire una scacchiera";
            }
        } else {
            std::cout << "La dimensione del reticolo non è multiplo della dimensione del blocco";
        }
        return 0;
    }


    int main() {
        using namespace std;
        unsigned seed = time(0);
        srand(seed);
        vector<int> lattice(N);
        const int numThread = setBlockSize(A);
        cout << numThread << endl;
        float energy = 0;
        int M = 0;
        initialize_lattice(lattice, energy, M);
        print_lattice(lattice);
        float T = 0.1;
        vector<float> results(1);
        results[0] = 1;

        std::vector<int> rand_vect;
        create_rand_vect(rand_vect);

        auto start = std::chrono::high_resolution_clock::now();  // Start timing before simulation
        while (T <= 0.2) {

            results.push_back(simulate(T, lattice, energy, M, rand_vect, numThread));

            print_lattice(lattice);
            T += 0.1;


            cout << abs(results.back()) << endl;
        }
        auto end = std::chrono::high_resolution_clock::now();  // End timing after simulation
        std::chrono::duration<double> elapsed = end - start;  // Calculate elapsed time
        cout << elapsed.count() << endl;


        return 0;

    }
