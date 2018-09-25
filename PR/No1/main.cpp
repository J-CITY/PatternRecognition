#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iomanip>

using namespace std;

int N = 0;
int ITERATION = 0;
int SIZE = 0;

class No1 {
    vector<int> boxes;
    double n = 0;
    int iter = 0;
    int Good = 0;
    double aRes = 0;
public:
    No1(int _n, int _iter) {
        n = _n;
        iter = _iter;
        boxes.resize(n);
    }

    void Run() {
        Good = 0;
        aRes = AnalyticalSolution();
        for (auto i = 0; i < iter; ++i) {
            Event();
        }
        Print();
    }

    double AnalyticalSolution() {
        double a = fact(n)/fact(1)/fact(n-1);
        double sum = 0;

        for (double i = 0; i < n-1; ++i) {
            sum += pow(-1, i) * fact(n-1)/fact(i)/fact(n-1-i) * pow(1-(1+i)/n, n);
        }
        return a*sum;
    }

    void Print() {
        double res = double(Good)/iter;
        cout << std::setw(12) << aRes << "|" << std::setw(12) << res << "|" << std::setw(12) << res-aRes << "\n";
    }

    void Event() {
        Clear();
        for (auto i = 0; i < n; ++i) {
            boxes[std::rand() % int(n)]++;
        }
        if (test()) {
            Good++;
        }
    }
private:
    long double fact(int N) {
        if(N < 0)
            return 0;
        if (N == 0)
            return 1;
        else
            return N * fact(N - 1);
    }
    void Clear() {
        for (auto i = 0; i < boxes.size(); ++i) {
            boxes[i] = 0;
        }
    }
    bool test() {
        int t = 0;
        for (auto i = 0; i < boxes.size(); ++i) {
            if (boxes[i] == 0) {
                t++;
            }
        }
        return t == 1 ? true : false;
    }
};

int main() {
    std::srand(unsigned(std::time(0)));

    cout << "Enter N, ITERATION, SIZE.\n";
    cin >> N >> ITERATION >> SIZE;

    cout << std::setw(12) << "Analytical" << "|" << std::setw(12) << "Program" << "|" << std::setw(12) << "R\n";
    No1 no1(N, ITERATION);
    for (auto i = 0; i < SIZE; ++i) {
        no1.Run();
    }
    return 0;
}
