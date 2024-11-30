#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <random>
#include <chrono>
#include <ctime>
#include <map>
#include <numeric>
using namespace std;

#include <omp.h>

//#define p1
//#define p2
#define p3
//#define p4




#define N 100000
const double EPSILON = 1e-6; // Ошибка

std::random_device rd;


//int generateRandomNumber() {
//	std::mt19937 mt(rd() | omp_get_thread_num());
//	std::uniform_int_distribution<> dist(1, 6);
//	return dist(mt); 
//}

#ifdef p1

int simulate(const int n) {
	std::mt19937 mt(rd() | omp_get_thread_num());
	std::uniform_int_distribution<> dist(0, 2 * n);
	int sum = 0;

	int* queue = new int[n * 2];
	for (int i = 0; i < n; i++) {
		queue[i] = 1;
		queue[n + i] = -1;
	}

	for (int i = 0; i < 2 * n; i++) {
		int j = dist(mt) % (2 * n);
		int temp = queue[j];
		queue[i] = queue[j];
		queue[j] = temp;
	}

	for (int i = 0; i < 2 * n; i++) {
		sum += queue[i];
		if (sum <= 0) { return 0; }
	}

	return 1;

}

int main()
{
	int n;
	cout << "Enter the number of consumers: ";
	cin >> n;


	int succesful = 0;

#pragma omp parallel num_threads(1)
	{
		cout << omp_get_num_threads();

#pragma omp for reduction(+:succesful)
		for (int i = 0; i < N; i++) {
			//cout << "Iteration " << i << endl;
			succesful += simulate(n);
		}

	}
	cout << "Consumers = " << n << " succesful = " << succesful << " N Simulation = " << N << endl;
	cout << "Probability = " << ((double)succesful) / (N) << " Probablity by theory = " << (double)(1.0 / (n + 1.0)) << " Error? = " << (((double)succesful / N - (double)(1.0 / (n + 1.0))) < EPSILON);
}
#endif // practice1

#ifdef p2

int generateRandomNumber() {
	std::mt19937 mt(rd() | omp_get_thread_num());
	std::uniform_int_distribution<> dist(1, 6);
	return dist(mt);
}



int main() {
	int countA = 0;
	int countB = 0;

#pragma omp parallel num_threads(1)
	{
		int local_countA = 0;
		int local_countB = 0;

#pragma omp for
		for (int i = 0; i < N; i++) {
			int die1 = generateRandomNumber();
			int die2 = generateRandomNumber();

			int sum = die1 + die2;
			if (sum % 2 == 0) {
				local_countB++;
				if (sum == 8) {
					local_countA++;
				}
			}
		}

#pragma omp atomic
		countA += local_countA;
#pragma omp atomic 
		countB += local_countB;
	}

	cout << static_cast<double> (countA) * std::pow(countB, -1) << endl;
}

#endif //p2

#ifdef p3

int main() {


	const int desk_size = 36;
	const int numOfAces = 4;

	int countA = 0, countB = 0, countA_given_B = 0;

#pragma omp parallel
	{
		std::mt19937 mt(rd() | omp_get_thread_num());

		int local_countA = 0;
		int local_countA_given_B = 0;
		int local_countB = 0;

		

#pragma omp for
		for (int i = 0; i < N; i++) {
			std::vector<int> desk(desk_size);
			std::fill(desk.begin(), desk.begin() + numOfAces, 1); // Ace
			std::fill(desk.begin() + numOfAces, desk.end(), 0); // Simple
			std::shuffle(desk.begin(), desk.end(), mt);

			if (desk[1] == 1) {
				local_countA++;
			}
			if (desk[0] == 1) {
				local_countB++;
				std::cout << desk.size() << std::endl;
				desk.erase(desk.begin());
				std::shuffle(desk.begin(), desk.end(), mt);
				std::cout << desk.size() << std::endl;
				if (desk[0] == 1) {
					local_countA_given_B++;
				}
			}
		}

#pragma omp atomic
		countA += local_countA;

#pragma omp atomic
		countB += local_countB;
#pragma omp atomic
		countA_given_B += local_countA_given_B;
	}

	double prob = static_cast<double> (countA) * std::pow(N, -1);
	double probGiv = static_cast<double> (countA_given_B) * std::pow(countB, -1);

	std::cout << "P(A) = " << prob << " P(A|B) = " << probGiv << " P(theory) = " << 3 * std::pow(35,-1) << std::endl;
}

#endif

#ifdef p4

// g++ -O3 -fopt-info-vec-all -march=native -openmp:experimental -o vectorize OpenMP-lab1.cpp

// SIMD
int main() {
	std::vector<float> a = { 2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0 };
	std::vector<float> b = { 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0 };
	//float result = std::inner_product(a.begin(), a.end(), b.begin(), 0.0f);
	//std::cout << "Результат (std::inner_product): " << result << std::endl;
	float result = 0.0f;
#pragma omp simd
	for (size_t i = 0; i < a.size(); ++i) {
		result += a[i] * b[i];
	}

	std::cout << "Результат (цикл): " << result << std::endl;
}

#endif