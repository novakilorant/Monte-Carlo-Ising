#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;


// Measure the Monte Carlo correlation time τ from the exponential decay of the average magnetization starting from the ferromagnetic “all spins up” configuration.

double tau() {}

// Modify the external field h in the range h ∈[0,10] and plot the M(h) curve.

double M() {}

// Measure the susceptibility using χ = ∂M/∂h|h=0. [Hint: Use small fields.]

double X() {}

// Verify the susceptibility from the fluctuations of the magnetization, χ = (〈M2〉−〈M〉2)/(kBT).

double X_M(vector<double> M, double k_BT) {
    double M_avg = 0.0;
    double M2_avg = 0.0;
    for (int m=0;m<M.size();m++){
        M_avg += m;
        M2_avg += m*m;
    }
    return 
}
