#include "gpe/solver.hpp"
#include "numerics.hpp"
#include "snapshot_archive.hpp"

namespace {
constexpr int OUTPUT_PRECISION = std::numeric_limits<double>::max_digits10;
constexpr double GRID_TOLERANCE = 1e-9;
}

gpe::SimulationResult gpe::GPE_solver(const SimulationParameters& parameters) {
    const double Length=parameters.length, delta_x=parameters.dx;
    const double delta_t=parameters.dt, t_init=parameters.t_initial, t_fin=parameters.t_final;
    const double m_tilde=parameters.mass_ratio, G_tilde=parameters.G;
    const double gamma=parameters.gamma, err_iter=parameters.picard_tolerance;
    const int iterations=parameters.max_picard_iterations, saving=parameters.config_save_every;
    const int velocity_saving=parameters.velocity_save_every;
    for (double value : {Length, delta_x, delta_t, t_init, t_fin, m_tilde, G_tilde,
                         gamma, err_iter, parameters.V_0, parameters.force})
        if (!isfinite(value)) throw runtime_error("Parameters must be finite.");
    const double points=Length/delta_x, steps=(t_fin-t_init)/delta_t;
    if (!(Length>0 && delta_x>0 && delta_t>0 && t_init>=0 && t_fin>t_init && gamma>0 && m_tilde>=0 && err_iter>0)
        || !isfinite(points) || points<3 || points>numeric_limits<int>::max()
        || abs(points-round(points))>GRID_TOLERANCE
        || !isfinite(steps) || steps<1 || steps>numeric_limits<int>::max()
        || abs(steps-round(steps))>GRID_TOLERANCE
        || iterations<1
        || saving<1
        || velocity_saving<1)
        throw runtime_error("Invalid grid, time range, or iteration settings; length/dx and (end-start)/dt must be integers.");
    const int space_points=static_cast<int>(round(points)), total_steps=static_cast<int>(round(steps));
    const int K_fin=static_cast<int>(iterations), save_every=static_cast<int>(saving);
    const int velocity_every=static_cast<int>(velocity_saving);
    Complex V_0(parameters.V_0,0), V=parameters.initial_velocity.value_or(V_0);
    VectorXcd Phi_n=VectorXcd::Ones(space_points), Delta_matrix(space_points);
    if (t_init>0 && (!parameters.initial_field.size() || !parameters.initial_velocity))
        throw runtime_error("Continuation requires initial_field and initial_velocity.");
    if (parameters.initial_field.size()) {
        if (parameters.initial_field.size()!=space_points || !parameters.initial_field.allFinite())
            throw runtime_error("Initial field must be finite and match the grid.");
        Phi_n=parameters.initial_field;
    }
    if (!isfinite(real(V)) || !isfinite(imag(V)))
        throw runtime_error("Initial velocity must be finite.");
    const filesystem::path output=parameters.output_directory;
    if (!filesystem::create_directory(output)) throw runtime_error("Output directory already exists.");
    ofstream velocity(output/"velocity.dat", ios::app), norm(output/"norm.dat"), convergence(output/"convergence.dat"), metadata(output/"parameters.txt");
    for (auto* stream : {&velocity,&norm,&convergence,&metadata}) {
        stream->exceptions(ios::badbit|ios::failbit);
        *stream<<setprecision(OUTPUT_PRECISION);
    }
    metadata<<"length "<<Length<<"\ndx "<<delta_x<<"\ndt "<<delta_t<<"\nmass_ratio "<<m_tilde
              <<"\ncoupling "<<G_tilde<<"\ngamma "<<gamma<<"\nvelocity_0 "<<real(V_0)
              <<"\nforce "<<parameters.force
              <<"\nstart "<<t_init<<"\nend "<<t_fin<<"\ntolerance "<<err_iter<<"\niterations "<<K_fin
              <<"\nsave_every "<<save_every<<"\nvelocity_every "<<velocity_every<<"\nrestart_velocity "<<real(V)<<"\nrestart_velocity_imag "<<imag(V)<<"\n";
    metadata.flush();
    velocity<<"# time real_velocity imaginary_velocity\n";
    norm<<"# time euclidean_norm dx_sum_abs_phi_squared\n";
    convergence<<"# time iterations iterate_difference_norm\n";
    auto record_velocity = [&](double time) {
        velocity<<time<<" "<<real(V)<<" "<<imag(V)<<"\n";
        velocity.flush();
    };
    SnapshotArchive archive(output/"configurations.tar.gz");
    auto record = [&](double time, int step) {
        norm<<time<<" "<<Phi_n.norm()<<" "<<delta_x*Phi_n.squaredNorm()<<"\n";
        const auto snapshot = output/("phi_"+to_string(step)+".dat");
        const auto temporary = filesystem::path(snapshot.string()+".tmp");
        ostringstream phi;
        phi.exceptions(ios::badbit|ios::failbit);
        phi<<setprecision(OUTPUT_PRECISION);
        for (int i=0;i<space_points;++i)
            phi<<i*delta_x-Length/2.0<<" "<<real(Phi_n(i))<<" "<<imag(Phi_n(i))<<"\n";
        archive.add(snapshot.filename().string(), phi.str());
        if (parameters.write_plain_snapshots) {
            ofstream plain(temporary);
            plain.exceptions(ios::badbit|ios::failbit);
            plain << phi.str();
            plain.close();
            filesystem::rename(temporary, snapshot);
        }
        velocity.flush();
        norm.flush();
        convergence.flush();
    };
    Delta_Function(Delta_matrix,space_points);
    Delta_matrix *= G_tilde/delta_x;
    Complex mult1=(0.25*I*delta_t*(1.0+m_tilde))/(delta_x*delta_x), mult2=0.5*delta_t*I, mult3;
    Complex A1,B1,C1;
    VectorXcd Phi_tilde(space_points),Phi_hat(space_points),Phi_k(space_points),Phi_n1(space_points),Absolute_Phi_tilde_adv(space_points),Absolute_Phi_fixed(space_points),Error_iter(space_points),Error_config(space_points),Temp(space_points);
    record_velocity(t_init);
    record(t_init,0);
    for (long long count=1;count<=total_steps;++count) {
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
         if (k > K_fin
             || !Phi_hat.allFinite() || !isfinite(Error_iter.norm()))
             {
                 ostringstream message;
                 message << setprecision(OUTPUT_PRECISION)
                         << "Nonlinear iteration failed at time " << t_init+count*delta_t
                         << " after " << min(k,K_fin) << " iterations; residual "
                         << Error_iter.norm() << ", tolerance " << err_iter;
                 throw runtime_error(message.str());
             }
         Phi_n1 = Phi_hat;//Phi_{n+1} obtained
         Error_config = Phi_n1 - Phi_n;//measuring difference between the two configs in two consecutive TIME STEPS (NOT ITERATION STEPS)
         //new velocity calculation and file
         V = Velocity_calculation(Phi_n,Phi_n1,space_points,
             V_0 + m_tilde*parameters.force*(t_init+count*delta_t),m_tilde,gamma);//Velocity writing
         
        Phi_n=Phi_n1;
        convergence<<t+delta_t<<" "<<min(k,K_fin)<<" "<<Error_iter.norm()<<"\n";
        convergence.flush();
        if (count%velocity_every==0 || count==total_steps) record_velocity(t_init+count*delta_t);
        if (count%save_every==0 || count==total_steps) record(t_init+count*delta_t,count);
    }
    archive.close();
    velocity.close();
    norm.close();
    convergence.close();
    return {output/"configurations.tar.gz", output/"velocity.dat", std::move(Phi_n), V, t_fin};
}
