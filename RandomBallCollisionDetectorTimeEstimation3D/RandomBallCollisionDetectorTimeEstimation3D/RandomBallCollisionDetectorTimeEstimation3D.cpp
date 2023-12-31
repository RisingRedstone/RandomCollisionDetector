// RandomBallCollisionDetectorTimeEstimation3D.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <math.h>
#include<vector>
#include <algorithm>
#include<iterator>
#include<chrono>
#include<fstream>



using namespace std;

typedef long double PrecType;

template<class T> void listify(T list, int sizes) {
    if (sizes) {
        cout << "[";
        for (int i = 0; i < sizes - 1; i++) {
            cout << list[i] << ", ";
        }
        cout << list[sizes - 1];
        cout << "]";
    }
}

template<class T> ostream& operator<<(ostream& out, vector<T> V) {
    listify(V, V.size());
    return out;
}

inline PrecType randomNum() {
    return static_cast <PrecType> (rand()) / static_cast <PrecType> (RAND_MAX);
}

inline int randIntRange(int min, int max) {
    return min + (rand() * (int)(max - min) / RAND_MAX);
}

class Ball {
private:
    PrecType x, y, z;
    PrecType radius;
public:
    Ball(PrecType x, PrecType y, PrecType z, PrecType radius) : x(x), y(y), z(z), radius(radius) {}
    Ball(PrecType radius) : radius(radius) {
        x = randomNum();
        y = randomNum();
        z = randomNum();
    }
    inline PrecType dist(Ball b2) {
        return ((x - b2.x) * (x - b2.x)) + ((y - b2.y) * (y - b2.y)) + ((z - b2.z) * (z - b2.z));
    }
    inline bool collide(Ball b2) {
        return ((x - b2.x) * (x - b2.x)) + ((y - b2.y) * (y - b2.y)) + ((z - b2.z) * (z - b2.z)) < ((b2.radius + radius) * (b2.radius + radius));
    }
    friend class Grid;
    friend ostream& operator<<(ostream& out, const Ball& b);
};

ostream& operator<<(ostream& out, const Ball& b) {
    cout << "{X: " << b.x << ", Y: " << b.y << ", Z: " << b.z << ", R: " << b.radius << " }";
    return out;
}

class Grid {
private:
    int N;
    int K;
    PrecType radius;
    vector<Ball>* grid;
    Ball** memoryAccumulate;
public:
    inline int coord(int x, int y, int z) {
        return x + (N * y) + (N * N * z);
    }
    Grid(int N, int K, PrecType radius) : N(N), K(K), radius(radius) {
        grid = new vector<Ball>[N * N * N];
        for (int i = 0; i < K; i++) {
            Ball b(radius);
            grid[coord(int((b.x * N) - (1e-5)), int((b.y * N) - (1e-5)), int((b.z * N) - (1e-5)))].push_back(b);
        }
        memoryAccumulate = new Ball * [K];
    }
    int calculateCols() {
        int counts = 0;
        int accumulated = 0;
        for (int i = 0; i < N - 1; i++) {
            for (int j = 0; j < N - 1; j++) {
                for (int k = 0; k < N - 1; k++) {

                    for (int x = 0; x < 2; x++) {
                        for (int y = 0; y < 2; y++) {
                            for (int z = 0; z < 2; z++) {
                                int iter = coord(i + x, j + y, k + z);
                                for (int n = 0; n < grid[iter].size(); n++) {
                                    memoryAccumulate[accumulated] = &grid[iter][n];
                                    accumulated++;
                                }
                            }
                        }
                    }

                    if (accumulated == 0) continue;
                    for (int l = 0; l < accumulated - 1; l++) {
                        for (int m = l + 1; m < accumulated; m++) {
                            counts += memoryAccumulate[l]->collide(*memoryAccumulate[m]);
                        }
                    }
                    accumulated = 0;
                }
            }
        }
        return counts;
    }
    friend ostream& operator<<(ostream& out, const Grid& G);
};

ostream& operator<<(ostream& out, const Grid& G) {
    int siz = (G.N * G.N * G.N);
    listify(G.grid, siz);
    return out;
}

class GridSamples {
private:
    vector<Grid> grids;
    int sampleNum;
public:
    GridSamples(int N, int K, PrecType radius, int sampleNum) : sampleNum(sampleNum) {
        for (int i = 0; i < sampleNum; i++) {
            grids.push_back(Grid(N, K, radius));
        }
    }
    double timeTestedRun() {
        double totTime = 0;
        chrono::time_point<chrono::steady_clock> start, end;
        for (int i = 0; i < grids.size(); i++) {
            start = chrono::steady_clock::now();
            grids[i].calculateCols();
            end = chrono::steady_clock::now();
            totTime += ((chrono::duration<double, milli>)(end - start)).count();
        }
        return totTime / (int)sampleNum;
    }
};

class TesterRun {
private:
    int Nlow, Nhigh, Nstep, Ngrad;
    int Klow, Khigh, Kstep, Kgrad;
    PrecType radius;
    int runs;
    int samplesPerRun;
    const char* fileName;
public:
    TesterRun(int Nlow, int Nhigh, int Nstep, int Klow, int Khigh, int Kstep, int runs, int samplesPerRun, PrecType radius, const char* fileName) :
        Nlow(Nlow), Nhigh(Nhigh), Nstep(Nstep), 
        Klow(Klow), Khigh(Khigh), Kstep(Kstep), 
        runs(runs), samplesPerRun(samplesPerRun), fileName(fileName), 
        radius(radius) {
        Ngrad = (Nhigh - Nlow) / Nstep;
        Kgrad = (Khigh - Klow) / Kstep;
    }

    void Run() {
        //Make this to run samples with random N and K and store the time values
        ofstream fout;
        int N, K;
        double timeStep;
        fout.open(fileName, ios::out);
        fout << "N,K,Time";
        for (int i = 0; i < runs; i++) {
            N = (randIntRange(0, Ngrad) * Nstep) + Nlow;
            K = (randIntRange(0, Kgrad) * Kstep) + Klow;
            timeStep = GridSamples(N, K, radius, samplesPerRun).timeTestedRun();
            fout << "\n" << N << "," << K << "," << timeStep;
        }
        fout.close();
    }

    int RunlowestForN(int K) {
        pair<int, double> mintimeStep = {-1, 1e10};
        for (int n = Nlow; n <= Nhigh; n++) {
            double timeStep = GridSamples(n, K, radius, samplesPerRun).timeTestedRun();
            cout << n << "\t" << timeStep <<endl;
            if (timeStep < mintimeStep.second) {
                mintimeStep = { n, timeStep };
            }
        }
        cout << mintimeStep.first << "\t" << mintimeStep.second;
        return mintimeStep.first;
    }

    void MapPoints() {
        ofstream fout;
        fout.open(fileName, ios::out);
        fout << "K,minimum N,Time";
        for (int k = Klow; k <= Khigh; k += Kstep) {
            cout << k << endl;
            for (int i = 0; i < runs; i++) {
                pair<int, double> mintimeStep = { -1, 1e10 };
                for (int n = Nlow; n <= Nhigh; n += Nstep) {
                    double timeStep = GridSamples(n, k, radius, samplesPerRun).timeTestedRun();
                    if (timeStep < mintimeStep.second) {
                        mintimeStep = { n, timeStep };
                    }
                }
                fout << "\n" << k << "," << mintimeStep.first << "," << mintimeStep.second;
            }
        }
        fout.close();
    }

    int NValue(int K) {
        double ratio = 0.14159534724294387;
        double Ktemp = (double)(K * (K-1));
        return int(pow(128 * ratio * Ktemp, ((double)1/(double)6)) + 0.5);
    }

    double PredictedTimeValue(int K) {
        double mult = 0.00020301698495638938 * 2.8746260481115466e-05;
        double Ktemp = (double)(K * (K - 1));
        return sqrt(128 * Ktemp * mult);
    }

    void MapPointsBasedOnPred() {
        ofstream fout;
        fout.open(fileName, ios::out);
        fout << "K,Time,Predicted Time";
        for (int k = Klow; k <= Khigh; k += Kstep) {
            int N = NValue(k);
            cout << k << "\t" << N << endl;
            double timeStep;
            double timePred = PredictedTimeValue(k);
            for (int i = 0; i < runs; i++) {
                timeStep = GridSamples(N, k, radius, samplesPerRun).timeTestedRun();
                fout << "\n" << k << "," << timeStep << "," << timePred;
            }
        }
        fout.close();
    }
};


// Made grid, now continue with creating many grids and then running the simulation :)

const int Nlow = 3, Nhigh = 15, Nstep = 1;
const int Klow = 20, Khigh = 1000, Kstep = 20;
const PrecType radius = 0.13;
const int runs = 2000;
const int samplesPerRun = 20;
int main()
{
    //TesterRun(Nlow, Nhigh, Nstep, Klow, Khigh, Kstep, runs, samplesPerRun, radius, "Data.csv").Run();
    TesterRun T(2, 10, 1, 
        40, 500, 10, 
        5, 250, 
        0.13, "MinTimeData.csv");
    T.MapPoints();
    //T.MapPointsBasedOnPred();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
