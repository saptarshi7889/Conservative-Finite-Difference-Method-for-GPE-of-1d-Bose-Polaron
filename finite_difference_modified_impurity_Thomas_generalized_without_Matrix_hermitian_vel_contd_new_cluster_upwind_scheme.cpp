//  Created by smajumdar on 17/04/2023.
//

//including Header files/libraries, to ensure all the functions and variables in the code works. Do not modify.
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <iostream>
#include<stdlib.h>
#include<math.h>
#include<cmath>
#include<random>
#include<iomanip>
#include"Eigen/Dense"
#include"Eigen/SVD"
#include<fstream>
#include<complex>
#include<stdio.h>
#include<bits/stdc++.h>

#define Pi 3.1415926535897

using namespace Eigen;
using namespace std;

typedef std::complex<double> Complex;
const Complex I(0, 1); //defines i=sqrt(-1).

//defining different functions required for the simulation. DO NOT ALTER ANYTHING IN THE FOLLOWING PART UNTIL AND UNLESS YOU WANT TO CHANGE SOMETHING IN THE CODE FUNDAMENTALLY

//This function is required to input previously simulated data as configuration file as initial configuration.
unordered_map<string,string> process_arguments(int argc, char** argv)
{
    std::string delimiter = "=";
    unordered_map<string,string> dict;
    for(int i = 1; i < argc; i++)
    {
        string str = string(argv[i]);
        size_t pos = str.find(delimiter);
        string name = str.substr(0,pos);
        string value = str.substr(pos+1);
        dict.insert(make_pair(name,value));
    }
    return dict;
}

//Thomas tridiagonal solver for solving the tridiagonal equation to obtain \phi^{n+1} on OPEN BOUNDARY CONDITION
void Thomas_tridiagonal_solver(Ref<VectorXcd> a_prime, complex<double> a_new, const Ref<const VectorXcd> b_n, complex<double> c_new, const Ref<const VectorXcd> d, int L)
{
    int i,j;
    VectorXcd c_prime(L-1), d_prime(L);
    //c_prime calculation
    for (i=0; i<L-1; ++i)
    {
        if(i==0)
        {
            c_prime(i) = c_new/b_n(i);
        }
        else
        {
            c_prime(i) = c_new/(b_n(i)-(a_new*c_prime(i-1)));
        }
    }
   //d_prime calculation
    for (i=0; i<L; ++i)
    {
        if(i==0)
        {
            d_prime(i) = d(i)/b_n(i);
        }
        else
        {
            d_prime(i) = (d(i) - (a_new*d_prime(i-1)))/(b_n(i) - (a_new*c_prime(i-1)));
        }
    }
  //matrix solution
    for(j=L-1; j>=0; j = j-1)
    {
    //cout<<"j"<<j<<endl;
        if(j==L-1)
            a_prime(j)=d_prime(j);
        else
            a_prime(j)=d_prime(j) - (c_prime(j)*a_prime(j+1));
    }
}

//Thomas tridiagonal solver for solving the tridiagonal equation to obtain \phi^{n+1} on periodic boundary condition, uses the tridiagonal solver function defined above
void Thomas_tridiagonal_solver_PBC(Ref<VectorXcd> a,  complex<double> a1, complex<double> b1, complex<double> c1, const Ref<const VectorXcd> d, int L)
{
    int i;
    complex<double> gamma_prime = -1.0*b1;
    VectorXcd u(L),v(L),b_new(L);
    //assign b_new
    for (i=0;i<L;++i)
    {
        if(i==0)
        {
            b_new(i)=b1-gamma_prime;
        }
        else if (i==L-1)
        {
            b_new(i) = b1 - ((c1*a1)/gamma_prime);
        }
        else
        {
            b_new(i)=b1;
        }
    }
    //assign u
    u = VectorXcd::Zero(L);
    u(0) = gamma_prime;
    u(L-1) = c1;
    //assign v
    v = VectorXcd::Zero(L);
    v(0) = 1.0+I*0.0;
    v(L-1) = a1/gamma_prime;
    //solve by Thomas algorithm
    VectorXcd y(L), q(L);
    Thomas_tridiagonal_solver(y,a1,b_new,c1,d,L);
    Thomas_tridiagonal_solver(q,a1,b_new,c1,u,L);
    //actual solution
    complex<double> x1,x2,x3;
    x1 = v.transpose()*y;
    x2 = v.transpose()*q;
    x2 = x2 + (1.0+0.0*I);
    x3 = x1/x2;
    a = y - x3*q;
}

//function to create the delta function potential at the impurity location
void Delta_Function(Ref<VectorXcd> a, int L) //function for the delta function matrix
{
a = VectorXcd::Constant(L,0.0+0.0*I);
a(L/2-1) = 1.0+0.0*I;
}

//takes a vector \phi(x) as input, and the output is |\phi(x)|^2, used in the GPE equation solution
void Non_Linear_Matrix(Ref<VectorXcd> a, const Ref<const VectorXcd> b, int L)
{
for (int i=0;i<L;++i)
     {
      a(i) = b(i)*conj(b(i));
     }
}

//function to calculate the velocity of the impurity at every time step
complex<double> Velocity_calculation(const Ref<const VectorXcd> b, const Ref<const VectorXcd> b_new, int L, complex<double> V_i, double m, double g)//to calculate velocity at every time step
{
complex<double> sum(0.0,0.0),sum_fin(0.0,0.0), sum1, sum2;
for (int i=0;i<L;++i)
  {
  if(i==L-1)
   {
       sum1=(conj(b(i))+conj(b_new(i)))*((b(0)-b(i))+(b_new(0)-b_new(i)));
       sum2=(b(i)+b_new(i))*((conj(b(0))-conj(b(i)))+(conj(b_new(0))-conj(b_new(i))));
       sum = sum+(sum1-sum2);
   }
  else
  {
      sum1=(conj(b(i))+conj(b_new(i)))*((b(i+1)-b(i))+(b_new(i+1)-b_new(i)));
      sum2=(b(i)+b_new(i))*((conj(b(i+1))-conj(b(i)))+(conj(b_new(i+1))-conj(b_new(i))));
      sum = sum+(sum1-sum2);
  }
  }
 sum_fin = V_i + (0.125*I*m/sqrt(g))*sum; // V(t) = V_0 + (i*m/sqrt(gamma))*sum
 return sum_fin;
}

// DO NOT ALTER ANYTHING IN THE ABOVE PART

//start of the main code

int main(int argc, char *argv[]) {
    
    //don't alter anything here, is required to input configuration file from previous simulation.
        string folder = "";
        if (argc >1)
        {
            cout<<"Trying to load the file \n";
            unordered_map<string,string> args = process_arguments(argc,argv);
                for (std::pair<string, string> element : args)
                {
                   if (element.first=="folder")
                   {
                       folder = element.second;
                   }
                }
        }

//define parameters

    double m_tilde=1.0/8.0, G_tilde=-1.5, gamma=0.1; //defining parameters, m_tilde=m/M, G_tilde, and gamma
    double delta_t=0.001, delta_x = (1.0/10000.0),t_init=142.5; //defining time and space steps, change as require. See that t_init is now defined HERE, this needs to match EXACTLY with the time stamp of the configuration input file.
    Complex mult1 = (0.25*I*delta_t*(1.0+m_tilde))/(delta_x*delta_x), mult2 = (0.5*delta_t*I), mult3; //defining the different multiplicative factors used in solving GPE equation. Check eqn. (4.6) and (4.7) in long PDF, mult1 and mult3 are the multiplicative factors of L1 and L2, respectively. mult2 is multiplicative factor of the last term, and mult3 will be defined later as it depends on impurity velocity V(t).
    complex<double> V_0(2.0,0.0),V;//defining the initial velocity V_0, V will be used to store impurity velocity at every time step.
    int Length=2000,i;  //Length denotes the system size
    int space_points = int(Length/delta_x);//calculating number of total spatial points
    double err_iter = 10e-9, err_config = 10e-5;
    cout<<"G: "<<G_tilde<<" delta_t: "<<delta_t<<" delta_x: "<<delta_x<<" Length: "<<Length<<" m_tilde: "<<m_tilde<<" V_0/v: "<<real(V_0)<<endl; //outputs the simulation parameters to the output file
    
    //time measurement, is required to stop the code in time
    clock_t start, end, current;//time measurement
    start = clock();

//Enter the main code
//initialization
    VectorXcd Phi_init(space_points),Delta_matrix(space_points),Phi_init_real(space_points),Phi_init_imag(space_points);
    Delta_Function(Delta_matrix,space_points);
    Delta_matrix = (G_tilde/delta_x)*Delta_matrix;//delta_function with G_tilde strength initialized
    
    //Initial configuration being initialized
    char filename_input_Phi[200];
    ifstream input_Phi;
    //taking Phi
    sprintf(filename_input_Phi,"Phi_solution_herm_G%gm%gV%gtime%gdt%gdx%gL%d_faster_finite_difference.dat",G_tilde,m_tilde,real(V_0),t_init,delta_t,delta_x,Length);// Initial configuration (data from last simulation opened)
    input_Phi.open((folder+filename_input_Phi).c_str());
    if ((!input_Phi)) {
        cout << "Cannot open file.\n";
        abort();
      }// If the previous filename is given wrong or is absent in the designated folder, the code will not run! double-check after submitting job!
    
    double temp_phi;

    for(i=0; i<3*space_points; ++i)
    {
      //Goes over line by line in the configuration file which has three columns. Column 1 is position, column 2 is Re(phi(x)), and column 3 is Im(phi(x)).
      input_Phi>>temp_phi;
      if(i%3==1)
      Phi_init_real(i/3) = temp_phi+0.0*I; //inputs the re(phi)
      if(i%3==2)
      Phi_init_imag(i/3) = temp_phi + 0.0*I; //inputs the im(phi)
    }
    Phi_init = Phi_init_real + I*Phi_init_imag; //phi_init is created here finally.
    //cout<<"Initial Phi is"<<endl<<Phi_init<<endl;
    
    //initialising the output files for data
    char filename_velocity[200], filename_solution[200],filename_norm[200], filename_convergence_check[200];
    ofstream myfile_velocity, myfile_solution,myfile_norm, myfile_convergence_check;
    int rep=9; //keeps track of the repetition number of the simulation with the same code.
    sprintf(filename_velocity,"impurity_velocity_herm_G%gm%gV%gdx%gdt%gL%d_upwind_faster_finite_difference_contd_rep%d_corrected.dat",G_tilde,m_tilde,real(V_0),delta_x,delta_t,Length,rep);
    sprintf(filename_norm,"Norm_evolution_herm_G%gm%gV%gdx%gdt%gL%d_upwind_faster_finite_difference_contd_rep%d_corrected.dat",G_tilde,m_tilde,real(V_0),delta_x,delta_t,Length,rep);
    sprintf(filename_convergence_check,"Convergence_check_upwind_G%gm%gV0%gdx%gdt%gL%d_faster_finite_difference_corrected.dat",G_tilde,m_tilde,real(V_0),delta_x,delta_t,Length);
    
    
    double t, t_fin =20000.0,t_safe=2760.0;//t is a moving variable to denote the time inside the final for loop. t_init is the initial time. t_fin is a very large final time that is never reached during the run time of the code. t_safe is the time IN MINUTES after which the simulation will be automatically stopped, maximum t_safe being 47 hours.
    int k, K_fin =12;//K_fin is the number of maximum iterations for iterative convergence of configurations, in general the configurations converge after 3-5 iterations.

    
    //time_loop starts
    VectorXcd Phi_tilde(space_points), Phi_hat(space_points),Phi_k(space_points),Phi_n(space_points), Phi_n1(space_points), Absolute_Phi_tilde_adv(space_points), Absolute_Phi_fixed(space_points),Error_iter(space_points),Error_config(space_points),Temp(space_points); //defining various matrices that will be used for GPE solution
    complex<double> A1,B1,C1; // complex constants relevant for the LHS of eq. (4.6) and (4.7), defined in more details later.
    Phi_n = Phi_init;//Phi at n=0 is the initial phi, in this case 1 for all sites.
    V=0.274999+I*0.0; //V at t=t_init, input manually here.
    
    //recording the norm of the system at t=0
    myfile_norm.open(filename_norm, ios::app);
    myfile_norm<<t_init<<" "<<Phi_n.norm()<<endl;
    myfile_norm.close();
    
    //recording the impurity velocity at t=0
    myfile_velocity.open(filename_velocity, ios::app);
    myfile_velocity<<t_init<<" "<<real(V)<<" "<<imag(V)<<endl;
    myfile_velocity.close();
    
    int count = 1; // defining a counter that will be used later in the time_loop
    
    for(t=t_init;t<=t_fin;t=t+delta_t)//main time_loop starts
     {
     mult3 = (0.5*delta_t*V)/delta_x; //multiplicative factor of L2 matrix is defined, updates every time step with the value of impurity velocity V
         
         //IMPLEMENTING THE UPWIND SCHEME HERE
         //A1,B1, and C1 are the multiplicative factors of LHS of eq. (4.7)
         if (real(V)<0) //When V^n is negative
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
             if(count % 50 ==0)
             {
                 cout<<"Now using upwind scheme "<<endl; //outputs every 50 time steps if the code has used upwind scheme
             }
                 
         }
         else if (real(V) >= 0) //When V is positive. A1,B1, and C1 are the multiplicative factors of LHS of eq. (4.6)
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
             } //Initializes the vector A_2 \times \Phi^n from the RHS in eqn. (4.6)
           if(count % 50 ==0)
            {
             cout<<"Now using normal scheme "<<endl; //outputs every 50 time steps (t+0.5) if the code has used upwind scheme
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
           Thomas_tridiagonal_solver_PBC(Phi_hat, A1, B1, C1, Phi_tilde, space_points); //Solves the equation A_1 Phi_{k+1}^{n+1} = Phi_tilde(i) to obtain Phi_{k+1}^{n+1}
           Error_iter = Phi_hat - Phi_k; //checks the differences between the configurations of two consecutive iteration step.
           
           //Checks convergence condition
           if(Error_iter.norm()<=err_iter)
               //when the norm of the error matrix is smaller than 10^{-9}, the iteration loop breaks and we store Phi^{n+1} = Phi_{k+1}^{n+1}. If not, goes in for another iteration loop to line 273.
           {
               if (count % 50 == 0) //outputs every 500 time steps at how many iterations the convergence condition is being met
               {
                   myfile_convergence_check.open(filename_convergence_check, ios::app);
                   myfile_convergence_check<<"Convergence at time "<<t<<"is at "<<k<<" th iteration with error "<<Error_iter.norm()<<endl;
                   myfile_convergence_check.close();
               }
               break;
           }
           else
           {
               Phi_k = Phi_hat;//Phi^{n+1}_k obtained
           }
       }
         Phi_n1 = Phi_hat;//Phi_{n+1} obtained
         Error_config = Phi_n1 - Phi_n;//measuring difference between the two configs in two consecutive TIME STEPS (NOT ITERATION STEPS)
         //new velocity calculation and file
         V = Velocity_calculation(Phi_n,Phi_n1,space_points,V_0,m_tilde,gamma);//Velocity writing
         
         //Density writing
         //1) checking if equilibrium has been reached. This condition is never matched during the entire runtime of the code because of the density wave never decaying to zero.
         if (Error_config.norm() <= err_config)
         {
             myfile_norm.open(filename_norm, ios::app);
             myfile_norm<<t+delta_t<<" "<<Phi_n1.norm()<<endl;//input Norm
             myfile_norm.close();
             myfile_velocity.open(filename_velocity,ios::app);
             myfile_velocity<<t+delta_t<<" "<<real(V)<<" "<<imag(V)<<endl;//input Velocity
             myfile_velocity.close();
             cout<<"final time for stationary configuration is : "<<t+delta_t<<endl;
             sprintf(filename_solution,"Phi_solution_G%gm%gV%gtime%gdt%gdx%gL%d_eqbm_finite_difference.dat",G_tilde,m_tilde,real(V_0),t+delta_t,delta_t,delta_x,Length);
             myfile_solution.open(filename_solution, fstream::out);
             for (int i=0;i<space_points;++i)
             {
                 myfile_solution<<(i*delta_x)-(Length/2.0)<<" "<<real(Phi_n1(i))<<" "<<imag(Phi_n1(i))<<endl;
             }
             myfile_solution.close();
             break;
         }
         
         //2) After certain time-interval, writing and saving impurity velocity and configuration file.
         else if(count % 100 == 0) //every 100 time step, so for (t+n*0.1), saves velocity and and norm. We usually save velocities in output file in 0.1 intervals.
         {
             myfile_norm.open(filename_norm, ios::app);
             myfile_norm<<t+delta_t<<" "<<Phi_n1.norm()<<endl;//output Norm
             myfile_norm.close();
             myfile_velocity.open(filename_velocity,ios::app);
             myfile_velocity<<t+delta_t<<" "<<real(V)<<" "<<imag(V)<<endl;//output Velocity
             myfile_velocity.close();
             if (count % 500 == 0) //every 250 time step, so for (t+n*0.5), saves configuration. For dt=0.01, we generally save configuration every 2.5 time step. For dt=0.0025, we save every 1 time step. for dt=0.001 and 0.0005, we save every 0.5 time steps.
             {
                 cout<<"time now is : "<<t+delta_t<<endl;
                 sprintf(filename_solution,"Phi_solution_herm_G%gm%gV%gtime%gdt%gdx%gL%d_faster_finite_difference.dat",G_tilde,m_tilde,real(V_0),t+delta_t,delta_t,delta_x,Length); //outputs Phi(x) with the TIME MENTIONED IN THE FILE NAME
                 myfile_solution.open(filename_solution, fstream::out);
                 for (int i=0;i<space_points;++i)
                 {
                     myfile_solution<<(i*delta_x)-(Length/2.0)<<" "<<real(Phi_n1(i))<<" "<<imag(Phi_n1(i))<<endl;
                 }
                 myfile_solution.close();
             }
             //IMPORTANT DUE TO CURRENT STRUCTURE OF THE CLUSTER.
             current = clock();//Checks the current time passed
             double time_passed = double(current-start)/double(CLOCKS_PER_SEC);
             if (double(time_passed/60.0) > double(t_safe))
                 break; // If time has passed, the code saves the last file and stops. IMPORTANT TO COPY ALL THE DATA FROM THE TMPDIR OF THE NODES. IF NOT DONE CORRECTLY, ALL FILES ARE LOST!
             Phi_n=Phi_n1;//new Phi_n is now Phi_{n+1}, used to be obtained next step.
             count++;
         }
        //3) continuing for the next time step if none of the previous time conditions is met.
        else
        {
            Phi_n=Phi_n1;//new Phi_n is now Phi_{n+1}, used to be obtained next step.
            count++;
        }
         
     }
    myfile_velocity.close();
    myfile_norm.close();
    end = clock();//time measurement stopped
    double time_taken = double(end-start)/double(CLOCKS_PER_SEC);
    cout<< "Time taken by the program is "<<time_taken/60.0<< " mins "<<endl;
    return 0;
}
