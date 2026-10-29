#include<iostream>
#include<tuple>
#include<vector>
#include<fstream>
#include<chrono>
#include<thread>
#include "capd/capdlib.h"
#include "Constants.h"
#include "cr3bp.h"
#include "MatrixTools.h"
#include "NonRigorousSolver.h"
#include "ComputeImageTask.h"

using namespace std;
using namespace capd;

// One non-rigorous (tau,e,x,y,z,dx,dy,dz) sample line from init_data.txt.
struct Record{
    long double tau;
    long double e;
    LDVector v;

    Record(const string& line) : v(6) {
        istringstream ss(line);
        if(!(ss >> tau >> e >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5]))
            throw runtime_error("Problem with parsing data");
    }
};

// A rigorous box for one boundary cell: the hull of two consecutive Records.
struct DataInstance{
    Interval tau;
    Interval E;
    IVector V;

    DataInstance(Interval tau, Interval E, IVector V) : tau(tau), E(E), V(V) {}
};

// init_data.txt holds 4 blocks of Records (the 4 walls of the (tau,E) rectangle's
// boundary, separated by blank lines.
// Each consecutive pair of Records within a block is hulled into one DataInstance,
// so the walls end up subdivided into many small boundary boxes.
vector<DataInstance> set_data(string file_name) {
    ifstream file(file_name);
    if(!file) throw runtime_error("File does not exist");

    vector<DataInstance> res;

    string line;

    for(int i = 0; i < 4; i++) {
        getline(file,line);
        if(line == "") break;
        Record left_record(line);

        while(getline(file,line)) {

            if(line == "") break;
            Record right_record(line);

            Interval tau(left_record.tau, right_record.tau);
            Interval E(left_record.e, right_record.e);
            IVector V = intervalHull(IVector(left_record.v),IVector(right_record.v));

            auto data = DataInstance(tau,E,V);
            res.push_back(data);

            left_record = right_record;
        }

    }
    return res;
}

// Builds boxes meant to cover the WHOLE (tau,E) rectangle,
// using only the left and right walls (3rd/4th blocks) of init_data.txt. 
// tau is the hull of the left wall's and right
// wall's (each constant) tau, i.e. the rectangle's full original tau-range. The right
// wall has E_div original boundary cells (E_div+1 Records), grouped into
// e_subdivision wider strips (group size = E_div/e_subdivision, so e_subdivision
// must divide e_div evenly).
vector<DataInstance> set_interior_data(string file_name, int e_subdivision) {
    ifstream file(file_name);
    if(!file) throw runtime_error("File does not exist");

    string line;
    // skip the first 2 blocks (bottom, top walls)
    for(int block = 0; block < 2; block++) {
        getline(file,line);
        while(getline(file,line)) if(line == "") break;
    }

    // left wall (3rd block): tau is constant across it, only need one sample
    getline(file,line);
    Record left_wall_sample(line);
    while(getline(file,line)) if(line == "") break;

    // right wall (4th block) in full
    vector<Record> right_wall;
    while(getline(file,line)) {
        if(line == "") break;
        right_wall.push_back(Record(line));
    }

    int E_div = right_wall.size() - 1; // number of original boundary cells
    if(E_div % e_subdivision != 0)
        throw runtime_error("e_subdivision must evenly divide the original E subdivision count");
    int group_size = E_div / e_subdivision; // original cells per new strip

    Interval tau(left_wall_sample.tau, right_wall[0].tau);

    vector<DataInstance> res;
    for(int g = 0; g < e_subdivision; g++) {
        const Record& first = right_wall[g * group_size];
        const Record& last = right_wall[(g + 1) * group_size];

        Interval E(first.e, last.e);
        IVector V = intervalHull(IVector(first.v), IVector(last.v));

        res.push_back(DataInstance(tau, E, V));
    }
    return res;
}

int main() {
    bool generate_file = 0;

    cout.precision(20);
    
    capd::threading::ThreadPool pool(threadsNo);
    init_cr3bp_pools();

    // half-widths and subdivision counts of the (tau,e) rectangle boundary around (TAU_TILDE,E_TILDE)
    NonRigorousSolver nonrig_solver(constants::DELTA_TAU,constants::DELTA_E,constants::TAU_DIV,constants::E_DIV);
    
    if(!ifstream("init_data.txt") || generate_file) {
        cout << "Generating file" << endl;
        nonrig_solver.save_init_data_to_file();
    } 

    // data for enclosing image of the boundary of the rectangle
    auto boundary_data = set_data("init_data.txt");
    
    // data for enclosing image of the interior of the rectangle
    int e_subdivision = (constants::E_DIV - 1) / 4;
    auto interior_data = set_interior_data("init_data.txt", e_subdivision);
    
    auto all_data = boundary_data;
    all_data.insert(all_data.end(), interior_data.begin(),interior_data.end());

    vector<ComputeImageTask*> tasks;
    // J_inv: non-rigorous isomorphism straightening the Poincare image near the
    // homoclinic point; J: its rigorous inverse, used by ComputeImageTask to
    // express the final image in the basis where the zero-finding argument applies
    IMatrix J_inv = nonrig_solver.get_isomorphism();
    IMatrix J = matrixAlgorithms::krawczykInverse(J_inv);

    auto start = chrono::high_resolution_clock::now();

    // output.txt: one line per boundary piece, in the order of init_data.txt, with the
    // enclosure of H([0,1], piece) as "x_left x_right y_left y_right" (the (dx,dz)
    // components after J). Interior strips are not saved. The file is written only
    // after all tasks finish.
    ofstream file("output.txt");
    file << std::scientific << std::setprecision(std::numeric_limits<long double>::max_digits10);
    
    // L's y-row/column are exactly zero (we work on the section Pi_y = {y=0}),
    // so L itself is singular -- invert the non-degenerate 5x5 submatrix, then reinsert.        
    IMatrix L = nonrig_solver.get_change_of_basis(0.);
    IMatrix L_temp = capd::matrixAlgorithms::krawczykInverse(matrix_erase_cord(L,1));
    IMatrix L_inv = matrix_add_cord(L_temp,1);
    
    int N = boundary_data.size();

    for(int i = 0; i < all_data.size(); i++) {

        Interval tau = all_data[i].tau;
        Interval E = all_data[i].E;

        IVector V = all_data[i].V;
        
        bool compute_boundary = (i < N);

        tasks.push_back(new ComputeImageTask(tau,E,V,J,L,L_inv,compute_boundary));
        
    }

    // Step 1 gives v(e) unique in the box V of each task. If [D F] is regular on the hull
    // of all these boxes, F(.,e) has at most one zero there, so all tasks give the same
    // v(e) and v is one C^1 curve on E. Checked here, before the tasks start, with the
    // solvers of thread 0.
    bool hull_ok = false;
    {
        IVector U = tasks[0]->get_V();
        for(auto* task: tasks) U = intervalHull(U, task->get_V());
        try {
            hull_ok = is_jacobian_regular(U,
                CR3BPPool::pool_Z->getSolver(0).pm, CR3BPPool::pool_E->getSolver(0));
        } catch(const std::exception& exc) {
            cout << exc.what() << endl;
        }
        if(!hull_ok) cout << "Jacobian is not regular on the hull of all boxes V" << endl;
    }
    cout << "Starting parallel computation with " << tasks.size() << " tasks:" << endl;
    for(auto* task: tasks) pool.process(task);

    // poll the shared counter and print a simple progress bar while tasks run
    int total = tasks.size();
    while(ComputeImageTask::completed_count < total) {
        int done = ComputeImageTask::completed_count;
        int barWidth = 40;
        int filled = total ? done * barWidth / total : barWidth;

        cout << "\r[" << string(filled,'#') << string(barWidth - filled,' ') << "] "
             << done << "/" << total << std::flush;

        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    cout << "\r[" << string(40,'#') << "] " << total << "/" << total << endl;

    for(auto* task: tasks) task->join();

    bool is_everything_ok = hull_ok;

    // (x,y) is an enclosure of a homotopy for (tau,E) box. If every instance of
    // computation gives status=COMPLETED, then the whole boundary surrounds zero
    // without touching it. Also it checks if G is defined in interior by evaluating
    // all the computations
    for(auto* task: tasks) {
        if(task->status != ComputeImageTask::COMPLETED) is_everything_ok = false;

        if(task->compute_boundary) {
            Interval x = task->result[0];
            Interval y = task->result[1];

            file << x.leftBound() << " " << x.rightBound() << " " 
                 << y.leftBound() << " " << y.rightBound() << endl;    
        }
        delete task;
    }

    auto stop = chrono::high_resolution_clock::now();
    auto total_seconds = chrono::duration_cast<chrono::seconds>(stop - start).count();
    cout << "Elapsed time: " << total_seconds / 60 << "min " 
         << total_seconds % 60 << "sec" << endl;

    if(is_everything_ok) {
        cout << "The existence of a homoclinic orbit proved successfully" << endl;
    }
    
    file.close();

    return 0;
}