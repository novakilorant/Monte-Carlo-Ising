#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>

// Pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


namespace py = pybind11;

using namespace std;

int n_rows, n_cols, time_steps;
string file_name;
enum initialState {UP, DOWN, RANDOM};

int seed = time(0);

struct MCIsing {
    int rows, cols;
    const double K, k_BT, h;
    vector<vector<int>> grid;

    MCIsing(int r, int c, initialState state = RANDOM, double K = 1, double k_BT = 5, double h = 0)
        : rows(r), cols(c),
          grid(r, vector<int>(c, 0)), K(K), k_BT(k_BT), h(h) {
        srand(seed);
        initialize(state);
    }

    void initialize(initialState state) {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j) {
                if (state == UP) grid[i][j] = 1;
                else if (state == DOWN) grid[i][j] = -1;
                else grid[i][j] = ((double)rand() / RAND_MAX < 0.5) ? 1 : -1;
            }
    }

    int getSpin(int i, int j) {
        if (i < 0) i = rows - 1;
        else if (i >= rows) i = 0;
        if (j < 0) j = cols - 1;
        else if (j >= cols) j = 0;
        return grid[i][j]; 
    }

    void step() {
        for (int i = 0; i < rows * cols; ++i) {
            int r = (int)rand() % rows;
            int c = (int)rand() % cols;
            Metropolis(r, c);
        }   
    }

    void display(ofstream& file) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (j == cols - 1) file << grid[i][j];
                else file << grid[i][j] << "\t";
            }
            file << '\n';
        }
    }

    void run_external_output(int time_steps, ofstream& file) {
        for (int t = 0; t < time_steps; ++t) {
            step();
        }
        display(file);
    }

    vector<vector<vector<int>>> run_numpy_output(int time_steps, ofstream& file) {
        vector<vector<vector<int>>> states;
        for (int t = 0; t < time_steps; ++t) {
            step();
            states.push_back(grid);
        }
        return states;
    }

    double H(const vector<vector<int>>& grid, int rows, int cols, double K, double k_BT, double h) {
        int adjacentSum = 0;
        int totalSum = 0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                adjacentSum += (getSpin(i+1, j) + getSpin(i-1, j) + getSpin(i, j+1) + getSpin(i, j-1)) * getSpin(i, j);
                totalSum += getSpin(i, j);
            }
        }
        return -K * adjacentSum + h * totalSum;
    }

    void Metropolis(int i, int j) {
        int relevantSum = getSpin(i+1, j) + getSpin(i-1, j) + getSpin(i, j+1) + getSpin(i, j-1);
        double deltaE = 2 * getSpin(i, j) * (K * relevantSum + h);
        if (deltaE <= 0) {
            grid[i][j] *= -1;
        } else {
            double prob = exp(-deltaE / k_BT);
            if (((double)rand() / RAND_MAX) < prob) {
                grid[i][j] *= -1;
            }
        }
    }
};



void getInput (int argc, char** argv) {

    if (argc != 7 ) {
        cout << " Usage 1: " << argv[0] << " time_steps n_rows n_cols initial_state K k_BT file_name" << endl;
    }

}

int main(int argc, char** argv) {

    getInput(argc, argv);
    ofstream file(file_name);

    MCIsing game(n_rows, n_cols);
    game.initialize();
    game.display(file);
    game.run(run_time-1, file);
    file.close();

    

    return 0;
}


PYBIND11_MODULE(example, m) {
    m.doc() = "C++ module for Monte Carlo Ising simulation"; // Optional module docstring
    m.def("MCIsing", [](int time_steps, int rows = 10, int cols = 10, initialState state = RANDOM, double K = 1, double k_BT = 5) {
        MCIsing sim(rows, cols, state, K, k_BT);
        std::ofstream file("output.txt");
        sim.run(time_steps, file);
        file.close();
    }, "A function that runs the MC Ising simulation");
