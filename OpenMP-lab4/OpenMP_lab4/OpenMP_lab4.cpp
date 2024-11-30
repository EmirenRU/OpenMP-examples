#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <random>
#include <chrono>
#include <ctime>
#include <map>
#include <cstring>
#include <numeric>
#include <iomanip>

using namespace std;

#include <omp.h>

#define p1
//#define p2

#define N 10000000
#define PI 3.14159265358979323846
const double EPSILON = 1e-6; // Ошибка

const int numberOfThreads = omp_get_max_threads() * 2;

std::random_device rd;
//std::mt19937 mt(rd());
//std::uniform_real_distribution<> dist(0, 6);

//double generateRandomNumber() { return dist(mt); }


#ifdef p1

double calculateMonteCarloVolume(double a, double b, int n, int threadsNum) {
    std::mt19937 mt(rd());
    std::uniform_real_distribution<> distX(0, b); // 0 < 𝑎 < 𝑏 -> till b?
    std::uniform_real_distribution<> distY(0, b); // same
    double vol = 0.0;
    int insideCounts = 0;
    #pragma omp parallel for reduction(+:vol, insideCounts) num_threads(threadsNum)
    for (int i = 0; i < n; i++) {

        double x = distX(mt);
        double y = distY(mt);

        if (x + y >= a && x + y <= b) {
            vol += sqrt(x * y);
            insideCounts++;
        }
    }

    // Sbigger - Ssmaller = S; (b^2 /2 - a^2/2)[a;b] + (a^2/2 - b^2/2)[-b;-a] = (b-a)(b+a) 
    double area = (b - a) * (b + a)  ;
    double average = vol * pow(insideCounts, -1);
    return average * area;
}


int main(){
    ofstream fout("results_1.csv");
    std::vector<double> times(numberOfThreads);
    std::vector<double> speedups(numberOfThreads);
    std::vector<double> efficiencies(numberOfThreads);

    double a = 1, b = 2;
    double appr = static_cast<double>(PI) * pow(12, -1) * (pow(b, 3) - pow(a, 3));

    for (int threads = 1; threads < numberOfThreads; threads++) {
        double t1 = omp_get_wtime();
        double result = calculateMonteCarloVolume(a, b, N, threads);
        double t2 = omp_get_wtime();
        double duration = t2 - t1;
        times[threads] = duration;

        if (threads == 1) {
            speedups[threads] = 1.0;
            efficiencies[threads] = 1.0;
        }
        else {
            speedups[threads] = times[1] / duration; // Using times[1] for speedup
            efficiencies[threads] = speedups[threads] / threads; // Efficiency
        }
        std::cout << "Threads: " << threads
            << ", Time: " << duration
            << ", Result: " << result
            << ", Approx: " << appr << std::endl;
        fout << threads << "," << duration << "," << speedups[threads] << "," << efficiencies[threads] << "\n";
    }
}



#endif

#ifdef p2

bool isInsideTheCircle(double x, double  y, double  x0, double  y0, double  r) {
    double dx = x - x0;
    double dy = y - y0;
    return (dx * dx + dy * dy) <= (r * r);
}

double calculateReuleauxTriangleArea(double r, int n, int threadNum) {
    
    int insideCounts = 0;
    #pragma omp parallel num_threads(threadNum)
    {
        std::mt19937 mt(rd() & omp_get_thread_num());
        std::uniform_real_distribution<> distX(-r, r);
        std::uniform_real_distribution<> distY(-r, r);
        #pragma omp for reduction(+: insideCounts) 
        for (int i = 0; i < n; i++) {
            double x = distX(mt);
            double y = distY(mt);

            //cout << " X = " << x << " Y = " << y << endl;
            /*if (isInsideTheCircle(x, y, r / sqrt(2), r / sqrt(2), r)) {
                std::cout << "Inside Circle 1: (" << x << ", " << y << ")\n";
            }
            if (isInsideTheCircle(x, y, -r / sqrt(2), r / sqrt(2), r)) {
                std::cout << "Inside Circle 2: (" << x << ", " << y << ")\n";
            }
            if (isInsideTheCircle(x, y, 0, -r, r)) {
                std::cout << "Inside Circle 3: (" << x << ", " << y << ")\n";
            }*/

            if (isInsideTheCircle(x, y, (1 / sqrt(2)), (1 / sqrt(2)), r) &&  // 45 degree
                isInsideTheCircle(x, y, (1 / sqrt(2)), (1 / sqrt(2)), r) && // 135 degree
                isInsideTheCircle(x, y, 0, -r, r)) // 270 degree
            {
                insideCounts++;
            }
        }
    }
    double squareArea = (2 * r) * (2 * r);
    cout << squareArea << " " << insideCounts << endl;

    return static_cast<double> (insideCounts) / n * squareArea;
}

int main() {
    ofstream fout("results_2.csv");
    std::vector<double> times(numberOfThreads);
    std::vector<double> speedups(numberOfThreads);
    std::vector<double> efficiencies(numberOfThreads);
    double r = 1.0;
    double approx = pow(2, -1) * (static_cast<double>(PI) - sqrt(3)) * pow(r, 2);
    for (int threads = 1; threads < numberOfThreads; threads++) {
        double t1 = omp_get_wtime();
        double result = calculateReuleauxTriangleArea(r, N, threads);
        double t2 = omp_get_wtime();
        double duration = t2 - t1;
        times[threads] = duration;

        if (threads == 1) {
            speedups[threads] = 1.0;
            efficiencies[threads] = 1.0;
        }
        else {
            speedups[threads] = times[1] / duration; // Using times[1] for speedup
            efficiencies[threads] = speedups[threads] / threads; // Efficiency
        }
        std::cout << "Threads: " << threads
            << ", Time: " << duration
            << ", Result: " << result
            << ", Approx: " << approx << std::endl;
        fout << threads << "," << duration << "," << speedups[threads] << "," << efficiencies[threads] << "\n";
    }

    fout.close();
}



#endif