#include "gpe/solver.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
using namespace std;
using namespace Eigen;
using Complex = std::complex<double>;

namespace {
constexpr double DEFAULT_LENGTH = 16.0;
constexpr double DEFAULT_DX = 0.125;
constexpr double DEFAULT_DT = 0.001;
constexpr double DEFAULT_END = 0.02;
constexpr double DEFAULT_MASS_RATIO = 0.125;
constexpr double DEFAULT_COUPLING = -1.5;
constexpr double DEFAULT_GAMMA = 0.1;
constexpr double DEFAULT_FORCE = 0.0;
constexpr double DEFAULT_VELOCITY = 2.0;
constexpr double DEFAULT_TOLERANCE = 1e-8;
constexpr int DEFAULT_ITERATIONS = 12;
constexpr int DEFAULT_SAVE_EVERY = 10;
constexpr double GRID_TOLERANCE = 1e-9;

void usage() {
    cout << "Usage: gpe [--length L] [--dx DX] [--dt DT] [--end T]\n"
         << "  [--mass-ratio m/M] [--coupling G] [--gamma GAMMA] [--velocity V0] [--force F]\n"
         << "  [--tolerance E] [--iterations K] [--save-every N] [--velocity-every N] [--output DIR]\n"
         << "  [--restart PHI_FILE --start T --restart-velocity V]\n"
         << "Restart requires the same physical and grid parameters, original V0,\n"
         << "and the velocity recorded at the checkpoint time. Output DIR must not exist.\n";
}
}

int main(int argc, char** argv) {
 try {
    map<string,string> args;
    const string allowed = " length dx dt end mass-ratio coupling gamma velocity tolerance iterations save-every velocity-every output restart start restart-velocity force ";
    for (int j=1; j<argc; ++j) {
        string key=argv[j];
        if (key=="--help") { usage(); return 0; }
        if (key.rfind("--",0)!=0 || allowed.find(" "+key.substr(2)+" ")==string::npos)
            throw runtime_error("Unknown option: "+key);
        if (++j==argc || args.count(key)) throw runtime_error("Missing or duplicate option: "+key);
        args[key]=argv[j];
    }
    auto number = [&](const string& key, double fallback) {
        if (!args.count(key)) return fallback;
        size_t used=0;
        double value=stod(args.at(key), &used);
        if (used!=args.at(key).size() || !isfinite(value)) throw runtime_error("Invalid number: "+key);
        return value;
    };
    const double Length=number("--length",DEFAULT_LENGTH), delta_x=number("--dx",DEFAULT_DX);
    const double delta_t=number("--dt",DEFAULT_DT), t_init=number("--start",0.0), t_fin=number("--end",DEFAULT_END);
    const double m_tilde=number("--mass-ratio",DEFAULT_MASS_RATIO), G_tilde=number("--coupling",DEFAULT_COUPLING);
    const double gamma=number("--gamma",DEFAULT_GAMMA), err_iter=number("--tolerance",DEFAULT_TOLERANCE);
    const double iterations=number("--iterations",DEFAULT_ITERATIONS), saving=number("--save-every",DEFAULT_SAVE_EVERY);
    const double velocity_saving=number("--velocity-every",saving);
    const double points=Length/delta_x, steps=(t_fin-t_init)/delta_t;
    if (!(Length>0 && delta_x>0 && delta_t>0 && t_init>=0 && t_fin>t_init && gamma>0 && m_tilde>=0 && err_iter>0)
        || !isfinite(points) || points<3 || points>numeric_limits<int>::max()
        || abs(points-round(points))>GRID_TOLERANCE
        || !isfinite(steps) || steps<1 || steps>numeric_limits<int>::max()
        || abs(steps-round(steps))>GRID_TOLERANCE
        || iterations<1 || iterations>numeric_limits<int>::max() || iterations!=floor(iterations)
        || saving<1 || saving>numeric_limits<int>::max() || saving!=floor(saving)
        || velocity_saving<1 || velocity_saving>numeric_limits<int>::max() || velocity_saving!=floor(velocity_saving))
        throw runtime_error("Invalid grid, time range, or iteration settings; length/dx and (end-start)/dt must be integers.");
    const int space_points=static_cast<int>(round(points));
    const int K_fin=static_cast<int>(iterations), save_every=static_cast<int>(saving);
    const int velocity_every=static_cast<int>(velocity_saving);
    Complex V_0(number("--velocity",DEFAULT_VELOCITY),0), V=V_0;
    VectorXcd Phi_n=VectorXcd::Ones(space_points), Delta_matrix(space_points);
    if (args.count("--restart")) {
        if (!args.count("--start") || !args.count("--restart-velocity"))
            throw runtime_error("Restart requires --start and --restart-velocity.");
        ifstream input(args.at("--restart"));
        for (int i=0;i<space_points;++i) {
            double x,re,im;
            if (!(input>>x>>re>>im) || !isfinite(x) || !isfinite(re) || !isfinite(im)
                || abs(x-(i*delta_x-Length/2.0))>GRID_TOLERANCE)
                throw runtime_error("Invalid checkpoint data or grid mismatch.");
            Phi_n(i)=Complex(re,im);
        }
        input>>ws;
        if (!input.eof()) throw runtime_error("Checkpoint has extra data.");
        V=Complex(number("--restart-velocity",0),0);
    } else if (args.count("--start") || args.count("--restart-velocity")) {
        throw runtime_error("--start and --restart-velocity require --restart.");
    }
    gpe::SimulationParameters settings;
    settings.length=Length; settings.dx=delta_x; settings.dt=delta_t;
    settings.t_initial=t_init; settings.t_final=t_fin;
    settings.mass_ratio=m_tilde; settings.G=G_tilde; settings.gamma=gamma;
    settings.picard_tolerance=err_iter; settings.max_picard_iterations=K_fin;
    settings.config_save_every=save_every; settings.velocity_save_every=velocity_every;
    settings.force=number("--force",DEFAULT_FORCE);
    settings.V_0=real(V_0); settings.initial_velocity=V; settings.initial_field=Phi_n;
    settings.output_directory=args.count("--output") ? args.at("--output") : "results";
    settings.write_plain_snapshots=true;
    gpe::GPE_solver(settings);
    if (args.count("--restart")) {
        ofstream metadata(settings.output_directory/"parameters.txt", ios::app);
        metadata.exceptions(ios::badbit|ios::failbit);
        metadata<<"restart "<<args.at("--restart")<<"\n";
        metadata.close();
    }
    cout<<"Results saved in "<<settings.output_directory<<"\n";
    return 0;
 } catch (const exception& error) {
    cerr<<"gpe: "<<error.what()<<"\n";
    return 1;
 }
}
