#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <random>

// Pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace std;

class MCIsing {
    public:
    const string file_name;
    string state;
    const int rows, cols;
    const double K, k_BT, h;
    vector<vector<int>> grid;
    vector<double> energy_record;
    vector<double> magnetization_record;
    int seed_random;
    mt19937 mt;
        MCIsing(int r, int c, string state = "RANDOM", double K = 1, double k_BT = 5, double h = 0, int seed = time(0), string file_name_ = "output.txt")
        : rows(r), cols(c),
          grid(r, vector<int>(c, 0)), K(K), k_BT(k_BT), h(h), state(state), file_name(file_name_), mt(seed) {
        initialize(state);
    }

    void initialize(string state = "RANDOM") {
        energy_record.clear();
        magnetization_record.clear();
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j) {
                if (state == "UP") grid[i][j] = 1;
                else if (state == "DOWN") grid[i][j] = -1;
                else grid[i][j] = (mt() % 2 == 0) ? 1 : -1;
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
            int r = mt() % rows;
            int c = mt() % cols;
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

    void run_external_output(int time_steps) {
        ofstream file(file_name);
        if (!file.is_open()) {
            cerr << "Error opening file: " << file_name << endl;
            return;
        }
        for (int t = 0; t < time_steps; ++t) {
            step();
            display(file);
            energy_record.push_back(H());
            magnetization_record.push_back(M());
        }
        file.close();
    }

    py::array_t<int> run_numpy_output(int time_steps) {
        auto result = py::array_t<int>({time_steps, rows, cols});
        auto states = result.mutable_unchecked<3>();
        for (int t = 0; t < time_steps; ++t) {
            step();
            // Copy the current state of the grid into the numpy array
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    states(t, i, j) = grid[i][j];
                }
            }
            energy_record.push_back(H());
            magnetization_record.push_back(M());
        }
        return result;
    }

    py::array_t<double> get_energy_record() {
        auto result = py::array_t<double>(energy_record.size());
        auto r = result.mutable_unchecked<1>();
        for (size_t i = 0; i < energy_record.size(); ++i) {
            r(i) = energy_record[i];
        }
        return result;
    }

    py::array_t<double> get_magnetization_record() {
        auto result = py::array_t<double>(magnetization_record.size());
        auto r = result.mutable_unchecked<1>();
        for (size_t i = 0; i < magnetization_record.size(); ++i) {
            r(i) = magnetization_record[i];
        }
        return result;
    }

    double H() {
        int adjacentSum = 0;
        int totalSum = 0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                adjacentSum += grid[i][j] * (getSpin(i+1, j) + getSpin(i-1, j) + getSpin(i, j+1) + getSpin(i, j-1));
                totalSum += grid[i][j];
            }
        }
        return -K * adjacentSum / 2.0 - h * totalSum;
    }

    double M() {
        int totalSum = 0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                totalSum += grid[i][j];
            }
        }
        return (double)totalSum / (rows * cols);
    }

    void Metropolis(int i, int j) {
        int relevantSum = getSpin(i+1, j) + getSpin(i-1, j) + getSpin(i, j+1) + getSpin(i, j-1);
        double deltaE = 2 * grid[i][j] * (K * relevantSum + h);
        uniform_real_distribution<double> dist(0.0, 1.0);
        if (deltaE <= 0 || dist(mt) < exp(-deltaE / k_BT)) {
            grid[i][j] *= -1;
        }
    }

    private:

};


PYBIND11_MODULE(simulator, m) {
    m.doc() = "Monte Carlo Ising Model Simulator using Metropolis Algorithm";
    py::class_<MCIsing>(m, "MCIsing")
        .def(py::init<int, int, string, double, double, double, int, string>(),
             py::arg("rows"),
             py::arg("cols"),
             py::arg("state") = "RANDOM",
             py::arg("K") = 1,
             py::arg("k_BT") = 5,
             py::arg("h") = 0,
             py::arg("seed") = time(0),
             py::arg("file_name") = "output.txt",
             "Initialize the MCIsing simulator with given parameters.\n\nParameters:\n----------\nrows : int\ncols : int\nstate : str\nK : float\nk_BT : float\nh : float\nseed : int\nfile_name : str\n")
        .def("step", &MCIsing::step)
        .def("initialize", &MCIsing::initialize, py::arg("state") = "RANDOM")
        .def("run_numpy_output", &MCIsing::run_numpy_output, py::arg("time_steps"))
        .def("run_external_output", &MCIsing::run_external_output, py::arg("time_steps"))
        .def("get_energy_record", &MCIsing::get_energy_record)
        .def("get_magnetization_record", &MCIsing::get_magnetization_record);
}
