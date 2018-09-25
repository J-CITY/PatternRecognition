#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iomanip>

using namespace std;

int ITERATION = 0;
int SIZE = 0;
int L = 0;

class No2 {
    double L = 0;
    int iter = 0;
    int Good = 0;
    double aRes = 0;
public:
    No2(double _L, int _iter) {
        iter = _iter;
        L = _L;
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
        double St = 0.5 * L/2 * L/2;
        double S = 0.5 * L * L;
        return St/S;
    }

    void Print() {
        double res = double(Good)/iter;
        cout << std::setw(12) << aRes << "|" << std::setw(12) << res << "|" << std::setw(12) << res-aRes << "\n";
    }

    void Event() {
        double a = 0;
        double b = 0;
        while(true) {
            a = drandom(0, L);
            b = drandom(0, L);
            if (a+b <= L) {
                break;
            }
        }
        if (a + b >= L/2 && a<=L/2 && b<=L/2) {
            Good++;
        }
    }
private:
    double drandom(double a, double b) {
        return a+(b-a)*std::rand()/RAND_MAX;
    }

};

int main() {
    std::srand(unsigned(std::time(0)));

    cout << "Enter L, ITERATION, SIZE.\n";
    cin >> L >> ITERATION >> SIZE;

    cout << std::setw(12) << "Analytical" << "|" << std::setw(12) << "Program" << "|" << std::setw(12) << "R\n";
    No2 no2(L, ITERATION);
    for (auto i = 0; i < SIZE; ++i) {
        no2.Run();
    }
    return 0;
}
