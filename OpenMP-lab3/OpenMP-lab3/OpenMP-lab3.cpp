#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <random>
#include <chrono>
#include <ctime>
#include <map>
#include <cstring>
#include <unordered_set>
#include <numeric>

using namespace std;

#include <omp.h>

//#define p1
//#define p2
//#define p3
//#define p4
#define p5 

#define N 1000000
const double REVERSED_N = pow(N, -1);
const double EPSILON = 1e-6; // Ошибка
const int numberOfThreads = omp_get_max_threads() * 2;

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_int_distribution<> *dist;
void initUniformRealDistribution(int a, int b) {
    dist = new uniform_int_distribution<>(a,b);
}

int generateRandomNumber() { return (*dist)(mt); }

int generateRandomNumber(std::mt19937 mt) { return (*dist)(mt); }



#ifdef p1

int main() {
    int a = 1, n = 6;
    initUniformRealDistribution(a, n);

    vector<int> distr(n,0);

    #pragma omp parallel num_threads(numberOfThreads)
    {
        std::mt19937 local_mt(omp_get_thread_num());
        vector<int> local_distr(n, 0);

        #pragma omp for 
        for (int i = 0; i <= N; i++) {
            int num = generateRandomNumber(local_mt);
            distr[num - 1]++;
        }

        #pragma omp critical
        {
            for (int i = 0; i < n; i++) {
                distr[i] += local_distr[i];
            }
        }
    }

    cout << "Distribution bentween 1 to " << n << " is ";
    for (int i = 0; i < n; i++) {
        cout << distr[i] << " ";
    }
    cout << endl << "Probability bentween 1 to " << n << " is ";
    for (int i = 0; i < n; i++) {
        cout << static_cast<double>(pow(N,-1)) * distr[i] << " ";
    }

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += pow(N, -1) * distr[i] ;
    }

    cout << endl << sum;

}

#endif //p1

#ifdef p2

int simulateMonteHall(bool switchChoice) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> dist(0,2);

    int prizeDoor = dist(mt);
    int playerChoice = dist(mt);

    if (!switchChoice) {
        return playerChoice == prizeDoor ? 1 : 0;
    }

    return playerChoice != prizeDoor ? 1 : 0;
}

int main() {
    initUniformRealDistribution(0, 2);

    int stayWins = 0;
    int switchWins = 0;


    #pragma omp parallel for num_threads(numberOfThreads) reduction(+:stayWins, switchWins)
    for (int i = 0; i < N; i++) {
        stayWins += simulateMonteHall(false);
        switchWins += simulateMonteHall(true);
    }

    double stayWinRate = static_cast<double>(stayWins) * REVERSED_N;
    double switchWinRate = static_cast<double>(switchWins) * REVERSED_N;

    cout << "Stay Win Rate = " << stayWinRate << endl;
    cout << "Switch Win Rate = " << switchWinRate << endl;

}
#endif

#ifdef p3

int simulateBirthday(int groupSize, int trials) {
    std::random_device rd;
    std::mt19937 mt(rd()&omp_get_thread_num());
    std::uniform_int_distribution<> dist(0, 364);

    int matches = 0;
    for (int trial = 0; trial < trials; trial++) {
        unordered_set<int> birthdays;
        for (int i = 0; i < groupSize; i++) {
            int birthday = dist(mt);
            if (birthdays.find(birthday) != birthdays.end()) {
                matches++;
                break;
            }
            birthdays.insert(birthday);
        }
    }
    return matches;
}


int main() {
    vector<double> freq(101, 0.0);
    ofstream fout("results.csv");

    #pragma omp parallel for num_threads(numberOfThreads)
    for (int n = 2; n < 100; n++) {
        cout << "Group: " << n << endl;
        int matches = simulateBirthday(n, N);
        freq[n] = static_cast<double>(matches) * REVERSED_N;
    }

    cout << "Group Size   |   Frequency" << endl;
    for (int i = 0; i < freq.size(); i++) {
        fout << i << "," << freq[i] << endl;
    }
    
}
#endif

#ifdef p4 

vector<double> simulateGamble(int both, int n1, int n2, double p) {
    std::mt19937 mt(rd() & omp_get_thread_num());
    std::uniform_real_distribution<double> dist(0, 1);
    //cout << dist(mt) << endl;

    int rounds = 0;

    int l1 = n1;
    int l2 = n2;

    while (l1 > 0 && l2 > 0) {
        rounds++;
        double randVal = dist(mt);
        if (randVal < p) {
            l1++;
            l2--;
        } else if (randVal >= p) {
            l1--;
            l2++;
        }
    }
    vector<double> v(3);
    v[0] = rounds;
    v[1] = l1;
    v[2] = l2;
    return v;
}

int main() {
    double p = 0.5;
    double q = 1 - p;
    int both = 100;
    int n1 = 50;
    int n2 = both - n1;


    double wins = 0, losses = 0;
    int totalRounds = 0;
    #pragma omp parallel for reduction(+:wins, losses, totalRounds)
    for (int i = 0; i < N; i++){
        vector<double> res = simulateGamble(both, n1, n2, p);
        totalRounds += (int) res[0];
        #pragma omp critical
        {
            cout << "Rounds: " << res[0] << " Wins: " << res[1] << " Losses: " << res[2] << endl;
        }
        if (res[1] == 0) {
            losses++;
        }
        else {
            wins++;
        }

    }
    std::cout << "Wins: " << wins << "\nLosses: " << losses << "\n";
    std::cout << "Average Rounds per game: " << static_cast<double>(totalRounds) * REVERSED_N << "\n";

}

#endif

#ifdef p5

double* simulateTokensGame(int x, int y, int z) {
    std::mt19937 mt(rd() ^ omp_get_thread_num());
    std::uniform_int_distribution<> dist(0, 2);
    cout << dist(mt) << endl;

    int rounds = 0;
    while (x >= 0 && y >= 0 && z >= 0) {
        rounds++;
        int winner = dist(mt);
        if (winner == 0) {
            x++;
            y--;
            z--;
        }
        else if (winner == 1) {
            x--;
            y++;
            z--;
        }
        else {
            x--;
            y--;
            z++;
        }
    }

    double* res = new double[4];
    res[0] = rounds;
    res[1] = x;
    res[2] = y;
    res[3] = z;
    return res;
}


int main() {
    int x, y, z;
    x = y = z = 10;
    int totalRounds = 0, w1 = 0, w2 = 0, w3 = 0;
    #pragma omp parallel for reduction(+:w1, w2, w3, totalRounds)
    for (int i = 0; i < N; i++) {
        double* res = simulateTokensGame(x,y,z);
        totalRounds += (int)res[0];
        if (res[1] > 0) {
            w1++;
        }
        else if (res[2] > 0){
            w2++;
        }
        else if (res[3] > 0) {
            w3++;
        }

        delete[] res;
    }
    std::cout << "Player 1: " << w1 << "\nPlayer 2: " << w2 << "\nPlayer 3: " << w3 << "\n";
    std::cout << "Average Rounds per game: " << static_cast<double>(totalRounds) * REVERSED_N << "\n";
}

#endif