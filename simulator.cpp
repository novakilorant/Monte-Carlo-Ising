#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <pybind11/pybind11.h>


namespace py = pybind11;

using namespace std;

int n_rows, n_cols, time_steps;
string file_name;
enum initialState {UP, DOWN, RANDOM};

int seed = time(0);

struct MCIsing {
    int rows, cols;
    vector<vector<int>> grid, newGrid;

    MCIsing(int r, int c)
        : rows(r), cols(c),
          grid(r, vector<int>(c, 0)), newGrid(r, vector<int>(c, 0)) {
        srand(seed);
    }

    void randomInitialize() {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                grid[i][j] = ((double)rand() / RAND_MAX < 0.5) ? 1 : -1;
    }

    int getCell(int i, int j) {
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

    void run(int generations, ofstream& file, double density) {
        for (int gen = 0; gen < generations; ++gen) {
            display(file);
            step(density);
        }
    }
};

void getInput (int argc, char** argv) {

    if (argc != 8 && argc != 6) {
        cout << " Usage 1: " << argv[0] << " n_rows n_cols n_rule boundary_condition density run_time file_name" << endl;
        cout << " n_rows: number of rows in the grid" << endl;
        cout << " n_cols: number of columns in the grid" << endl;
        cout << " n_rule: number of live neighbors to keep a cell alive" << endl;
        cout << " boundary_condition: OPEN, PERIODIC, ALIVE, RANDOM" << endl;
        cout << " density: initial density of live cells (0.0 to 1.0)" << endl;
        cout << " run_time: number of generations to run" << endl;
        cout << " file_name: name of the output file" << endl;
        cout << " This program simulates the Game of Life." << endl;
        cout << " Usage 2: " << argv[0] << " n_rows n_cols n_rule run_time file_name" << endl;
        cout << " n_rows: number of rows in the grid" << endl;
        cout << " n_cols: number of columns in the grid" << endl;
        cout << " n_rule: number of live neighbors to keep a cell alive" << endl;
        cout << " file_name: identifing name of the output file" << endl;
        cout << " This program simulates the Sand Dune model." << endl;
        exit(1);
    }

    if (argc == 6) {
        n_rows = atoi(argv[1]);
        n_cols = atoi(argv[2]);
        n_rule = atoi(argv[3]);
        run_time = atoi(argv[4]);
        file_name = argv[5];
        simulation_type = "SandDune";
        if (n_rule < 1) {
            cout << "n_rule must be greater than 0." << endl;
            exit(1);
        }
        return;
    }
    else if (argc == 8) {
        n_rows = atoi(argv[1]);
        n_cols = atoi(argv[2]);
        n_rule = atoi(argv[3]);
        string boundary_condition = argv[4];
        density = atof(argv[5]);
        if (boundary_condition == "OPEN") {
            boundary = OPEN;
        } else if (boundary_condition == "PERIODIC") {
            boundary = PERIODIC;
        } else if (boundary_condition == "ALIVE") {
            boundary = ALIVE;
        } else if (boundary_condition == "RANDOM") {
            boundary = RANDOM;
        } else {
            cout << "Invalid boundary condition. Use OPEN, PERIODIC, ALIVE, or RANDOM." << endl;
            exit(1);
        }
        run_time = atoi(argv[6]);
        file_name = argv[7];
        if (n_rule < 1) {
            cout << "n_rule must be greater than 0." << endl;
            exit(1);
        }
    }
}

int main(int argc, char** argv) {

    getInput(argc, argv);

    if (simulation_type == "GameOfLife") {
        ofstream file(file_name);

        MCIsing game(n_rows, n_cols);
        game.initialize();
        game.display(file);
        game.run(run_time-1, file);
        file.close();
    }
    

    return 0;
}


PYBIND11_MODULE(example, m) {
    m.doc() = "C++ module for Monte Carlo Ising simulation"; // Optional module docstring
    m.def("MCIsing", [](int time_steps, int rows = 10, int cols = 10, initialState state = RANDOM) {
        MCIsing sim(rows, cols);
        std::ofstream file("output.txt");
        sim.run(time_steps, file, state);
        file.close();
    }, "A function that runs the MC Ising simulation");
