#include "numerics.hpp"

namespace {
constexpr double DEFAULT_LENGTH = 16.0;
constexpr double DEFAULT_DX = 0.125;
constexpr double DEFAULT_DT = 0.001;
constexpr double DEFAULT_END = 0.02;
constexpr double DEFAULT_MASS_RATIO = 0.125;
constexpr double DEFAULT_COUPLING = -1.5;
constexpr double DEFAULT_GAMMA = 0.1;
constexpr double DEFAULT_VELOCITY = 2.0;
constexpr double DEFAULT_TOLERANCE = 1e-8;
constexpr int DEFAULT_ITERATIONS = 12;
constexpr int DEFAULT_SAVE_EVERY = 10;
constexpr int OUTPUT_PRECISION = std::numeric_limits<double>::max_digits10;
constexpr double GRID_TOLERANCE = 1e-9;

void usage() {
    cout << "Usage: gpe [--length L] [--dx DX] [--dt DT] [--end T]\n"
         << "  [--mass-ratio m/M] [--coupling G] [--gamma GAMMA] [--velocity V0]\n"
         << "  [--tolerance E] [--iterations K] [--save-every N] [--output DIR]\n"
         << "  [--restart PHI_FILE --start T --restart-velocity V]\n"
         << "Restart requires the same physical and grid parameters, original V0,\n"
         << "and the velocity recorded at the checkpoint time. Output DIR must not exist.\n";
}
}

int main(int argc, char** argv) {
 try {
    map<string,string> args;
    const string allowed = " length dx dt end mass-ratio coupling gamma velocity tolerance iterations save-every output restart start restart-velocity ";
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
    const double points=Length/delta_x, steps=(t_fin-t_init)/delta_t;
    if (!(Length>0 && delta_x>0 && delta_t>0 && t_init>=0 && t_fin>t_init && gamma>0 && m_tilde>=0 && err_iter>0)
        || !isfinite(points) || points<3 || points>numeric_limits<int>::max()
        || abs(points-round(points))>GRID_TOLERANCE
        || !isfinite(steps) || steps<1 || steps>numeric_limits<int>::max()
        || abs(steps-round(steps))>GRID_TOLERANCE
        || iterations<1 || iterations>numeric_limits<int>::max() || iterations!=floor(iterations)
        || saving<1 || saving>numeric_limits<int>::max() || saving!=floor(saving))
        throw runtime_error("Invalid grid, time range, or iteration settings; length/dx and (end-start)/dt must be integers.");
    const int space_points=static_cast<int>(round(points)), total_steps=static_cast<int>(round(steps));
    const int K_fin=static_cast<int>(iterations), save_every=static_cast<int>(saving);
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
    const filesystem::path output=args.count("--output") ? args.at("--output") : "results";
    if (!filesystem::create_directory(output)) throw runtime_error("Output directory already exists.");
    ofstream velocity(output/"velocity.dat"), norm(output/"norm.dat"), convergence(output/"convergence.dat"), parameters(output/"parameters.txt");
    for (auto* stream : {&velocity,&norm,&convergence,&parameters}) {
        stream->exceptions(ios::badbit|ios::failbit);
        *stream<<setprecision(OUTPUT_PRECISION);
    }
    parameters<<"length "<<Length<<"\ndx "<<delta_x<<"\ndt "<<delta_t<<"\nmass_ratio "<<m_tilde
              <<"\ncoupling "<<G_tilde<<"\ngamma "<<gamma<<"\nvelocity_0 "<<real(V_0)
              <<"\nstart "<<t_init<<"\nend "<<t_fin<<"\ntolerance "<<err_iter<<"\niterations "<<K_fin
              <<"\nsave_every "<<save_every<<"\nrestart_velocity "<<real(V)<<"\n";
    if (args.count("--restart")) parameters<<"restart "<<args.at("--restart")<<"\n";
    velocity<<"# time real_velocity imaginary_velocity\n";
    norm<<"# time euclidean_norm dx_sum_abs_phi_squared\n";
    convergence<<"# time iterations iterate_difference_norm\n";
    auto record = [&](double time, int step) {
        velocity<<time<<" "<<real(V)<<" "<<imag(V)<<"\n";
        norm<<time<<" "<<Phi_n.norm()<<" "<<delta_x*Phi_n.squaredNorm()<<"\n";
        ofstream phi(output/("phi_"+to_string(step)+".dat"));
        phi.exceptions(ios::badbit|ios::failbit);
        phi<<setprecision(OUTPUT_PRECISION);
        for (int i=0;i<space_points;++i)
            phi<<i*delta_x-Length/2.0<<" "<<real(Phi_n(i))<<" "<<imag(Phi_n(i))<<"\n";
    };
    Delta_Function(Delta_matrix,space_points);
    Delta_matrix *= G_tilde/delta_x;
    Complex mult1=(0.25*I*delta_t*(1.0+m_tilde))/(delta_x*delta_x), mult2=0.5*delta_t*I, mult3;
    Complex A1,B1,C1;
    VectorXcd Phi_tilde(space_points),Phi_hat(space_points),Phi_k(space_points),Phi_n1(space_points),Absolute_Phi_tilde_adv(space_points),Absolute_Phi_fixed(space_points),Error_iter(space_points),Error_config(space_points),Temp(space_points);
    record(t_init,0);
    for (int count=1;count<=total_steps;++count) {
        const double t=t_init+(count-1)*delta_t;
        int i,k;
     mult3 = (0.5*delta_t*V)/delta_x; //multiplicative factor of L2 matrix is defined, updates every time step with the value of impurity velocity V
     
         //IMPLEMENTING THE UPWIND SCHEME HERE
         //A1,B1, and C1 are the multiplicative factors of LHS of eq. (4.7)
         if (real(V)<0) //When V is negative
         {
             B1 = 1.0 + (2.0*mult1) - mult3; // for Phi_{k+1}^{n+1}(i,i), diagonal elements.
             A1= -mult1 + mult3; // for Phi_{k+1}^{n+1}(i,i+1), (diagonal+1) elements.
             C1 = -mult1; // for Phi_{k+1}^{n+1}(i+1,i), (diagonal-1) elements.
             for (i=0; i < space_points; ++i)
             {
                 if(i==0)
                 {
                     Temp(i) = mult1*Phi_n(i+1) + (mult1-mult3)*Phi_n(space_points-1) + (1.0-2.0*mult1+mult3)*Phi_n(i);
                 }
                 else if (i==space_points-1)
                 {
                     Temp(i) = mult1*Phi_n(0) + (mult1-mult3)*Phi_n(i-1) + (1.0-2.0*mult1+mult3)*Phi_n(i);
                 }
                 else
                 {
                     Temp(i) = mult1*Phi_n(i+1) + (mult1-mult3)*Phi_n(i-1) + (1.0-2.0*mult1+mult3)*Phi_n(i);
                 }
             } //Initializes the vector A_2 \times \Phi^n from the RHS in eqn. (4.7)
             if(count % 500 ==0)
             {
                 cout<<"Now using upwind scheme "<<endl; //outputs every 500 time steps if the code has used upwind scheme
             }
                 
         }
         else if (real(V) >= 0)//When V is positive. A1,B1, and C1 are the multiplicative factors of LHS of eq. (4.6)
         {
             B1 = 1.0 + (2.0*mult1) + mult3; // for Phi_{k+1}^{n+1}(i,i), diagonal elements.
             A1= -mult1; // for Phi_{k+1}^{n+1}(i,i+1), (diagonal+1) elements.
             C1 = -(mult1+mult3); // for Phi_{k+1}^{n+1}(i+1,i), (diagonal-1) elements.
             for (i=0; i < space_points; ++i)
             {
                 if(i==0)
                 {
                     Temp(i) = (mult1+mult3)*Phi_n(i+1) + mult1*Phi_n(space_points-1) + (1.0-2.0*mult1-mult3)*Phi_n(i);
                 }
                 else if (i==space_points-1)
                 {
                     Temp(i) = (mult1+mult3)*Phi_n(0) + mult1*Phi_n(i-1) + (1.0-2.0*mult1-mult3)*Phi_n(i);
                 }
                 else
                 {
                     Temp(i) = (mult1+mult3)*Phi_n(i+1) + mult1*Phi_n(i-1) + (1.0-2.0*mult1-mult3)*Phi_n(i);
                 }
             }//Initializes the vector A_2 \times \Phi^n from the RHS in eqn. (4.6)
           if(count % 500 ==0)
            {
             cout<<"Now using normal scheme "<<endl; //outputs every 500 time steps (t+0.5) if the code has used upwind scheme
            }
         }
         //A1 and A2*Phi_n assigned based on sign of V
         
         Phi_k = Phi_n;//before every picard iteration, set Phi_k to Phi_n
         Non_Linear_Matrix(Absolute_Phi_fixed,Phi_n,space_points);//|Phi_n|^2 matrix for every time step which is fixed during all the iterations, written in the notes as Abs(\Phi^n)^2.
         
         //Picard iteration for calculating Phi^{n+1} from Phi^n iteratively
      for(k=1;k<=K_fin;++k)
       {
           Non_Linear_Matrix(Absolute_Phi_tilde_adv,Phi_k,space_points); //Calculates |Phi_k^{n+1}|^2 at every convergence step from equation (4.6/4.7).
           
           //Calculates the RHS of equation (4.6/4.7)
           for (i=0; i < space_points; ++i)
           {
               Phi_tilde(i) = Temp(i) - mult2*(Delta_matrix(i) + 0.5*(Absolute_Phi_tilde_adv(i)+Absolute_Phi_fixed(i)))*(Phi_k(i)+Phi_n(i));
           }
           Thomas_tridiagonal_solver_PBC(Phi_hat, A1, B1, C1, Phi_tilde, space_points);//Solves the equation A_1 Phi_{k+1}^{n+1} = Phi_tilde(i) to obtain Phi_{k+1}^{n+1}
           Error_iter = Phi_hat - Phi_k; //checks the differences between the configurations of two consecutive iteration step.
           
           //Checks convergence condition
           if(Error_iter.norm()<=err_iter) //when the norm of the error matrix is smaller than 10^{-9}, the iteration loop breaks and we store Phi^{n+1} = Phi_{k+1}^{n+1}. If not, goes in for another iteration loop to line 273.
           {
               break;
           }
           else
           {
               Phi_k = Phi_hat;//Phi^{n+1}_k obtained
           }
       }
         if (k > K_fin || !Phi_hat.allFinite())
             throw runtime_error("Nonlinear iteration failed at time " + to_string(t));
         Phi_n1 = Phi_hat;//Phi_{n+1} obtained
         Error_config = Phi_n1 - Phi_n;//measuring difference between the two configs in two consecutive TIME STEPS (NOT ITERATION STEPS)
         //new velocity calculation and file
         V = Velocity_calculation(Phi_n,Phi_n1,space_points,V_0,m_tilde,gamma);//Velocity writing
         
        Phi_n=Phi_n1;
        convergence<<t+delta_t<<" "<<k<<" "<<Error_iter.norm()<<"\n";
        if (count%save_every==0 || count==total_steps) record(t_init+count*delta_t,count);
    }
    cout<<"Results saved in "<<output<<"\n";
    return 0;
 } catch (const exception& error) {
    cerr<<"gpe: "<<error.what()<<"\n";
    return 1;
 }
}
