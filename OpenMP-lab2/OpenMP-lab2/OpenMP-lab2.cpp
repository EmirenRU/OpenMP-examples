#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <random>
#include <chrono>
#include <ctime>
#include <cstring>
#include <omp.h>
#include <stdio.h>
using namespace std;


//#define p1
#define p2

#define N 10000000
const double EPSILON = 1e-6; // Ошибка

const int numberOfThreads = omp_get_max_threads() * 2;

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<> dist(0, 6);


double generateRandomNumber() {
    std::mt19937 mt(rd());
    std::uniform_real_distribution<> dist(0, 6);
    return dist(mt);

}


#ifdef p1

double omp_sum(vector<double>& v, int num_thr) {
    double sum = 0.0;
    size_t length = v.size();
    int threads = num_thr;

#pragma omp parallel for reduction(+:sum) num_threads(threads)
    for (int i = 0; i < length; i++)
    {
        sum += v[i];
    }
    return sum;
}

double omp_max(vector<double>& v, int num_thr) {
    double mV = v[0];
    size_t length = v.size();
    int threads = num_thr;


#pragma omp parallel for reduction(max:mV) num_threads(threads)
    for (int i = 0; i < length; i++) {
        mV = std::max(v[i], mV);
    }

    return mV;
}

double omp_min(vector<double>& v, int num_thr) {
    double mV = v[0];
    size_t length = v.size();
    int threads = num_thr;

#pragma omp parallel for reduction(min:mV) num_threads(threads)
    for (int i = 0; i < length; i++) {
        mV = std::min(v[i], mV);
    }

    return mV;
}

std::ofstream createFileForOutput(std::string name) {
    std::ofstream outFile(name);
    if (!outFile) {
        cerr << "Something went wrong with the file";
    }
    return outFile;
}


int main() {
    std::ofstream outFile;

    std::vector<double> array(N);


#pragma omp parallel for num_threads(omp_get_max_threads())
    for (int i = 0; i < array.size(); i++) {
        array[i] = generateRandomNumber();
    }

    int opt;
    cout << "What you want to do? \n1. Sum\n2. Max\n3. Min\nWrite a number (Example: 1): ";
    cin >> opt;

    double (*func)(vector<double>&, int);

    switch (opt) {
    case 1:
        outFile = createFileForOutput("result_sum.csv");
        cout << "Testing \"SUM\" method\n";
        func = &omp_sum;
        break;
    case 2:
        outFile = createFileForOutput("result_max.csv");
        cout << "Testing \"MAX\" method\n";
        func = &omp_max;
        break;
    case 3:
        outFile = createFileForOutput("result_min.csv");
        cout << "Testing \"MIN\" method\n";
        func = &omp_min;
        break;
    default:
        return 1;
    }
    const double test_num = 10;
    double avgTime = 0.0;

    outFile << "Threads,avgTime" << endl;
    for (int num = 0; num < numberOfThreads; num++) {
        double t1 = omp_get_wtime();
        double arr_res = (*func)(array, num);
        double t2 = omp_get_wtime();
        avgTime = (t2 - t1);
        cout << "Threads: " << num << " \n avgTime: " << avgTime * std::pow(2, -1) << endl;
        outFile << num + 1 << "," << avgTime << endl;
    }

}

#endif

#ifdef p2

double fi(double x) {
    return static_cast<double>(pow(x, 2) * 0.5);
}

double f(double x) {
    return x; // integral_1? res = x**2 / 2.0
}

double trapezoidal(double (*func)(double), double a, double b, int n, int threads_num) {
    double h = (b - a) / n;
    double sum = 0.0;
    double fa = 0, fb = 0, res = 0;

#pragma omp parallel for num_threads(threads_num) reduction(+:sum)
    for (int i = 0; i < n; i++) {
        sum += fi(a + i * h);
    }

    //#pragma omp atomic
    //fa += f(a);
    //#pragma omp atomic
    //fb += f(b);
    //#pragma omp atomic
    //res += h * 0.5 * (fa + fb) + sum * h;

    return h * 0.5 * (f(a) + f(b)) + sum * h;;
}

double get_analytic_solution(double a, double b) {
    return ((b * b) / 2) - ((a * a) / 2);
}



int main() {
    std::ofstream outFile("results_2.csv");
    if (!outFile) {
        cerr << "Something went wrong with the file";
        return 1;
    }
    double a = 0, b = 3;
    int iterations = N;
    double avgTime = 0.0;
    outFile << "Threads,avgTime" << endl;

    for (int numThread = 0; numThread < numberOfThreads; numThread++) {
        double t1 = omp_get_wtime();
        double res = trapezoidal((*f), a, b, iterations, numThread);
        double t2 = omp_get_wtime();
        avgTime = (t2 - t1);
        cout << "Threads: " << numThread << " \n avgTime: " << avgTime << endl;
        outFile << numThread << "," << avgTime << endl;
        cout << "Result is " << res << " Analytic solution is " << get_analytic_solution(a, b) << endl;
    }

}
#endif