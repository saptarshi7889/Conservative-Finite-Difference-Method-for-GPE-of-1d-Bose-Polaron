#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

// Numerical routines preserved from the supplied fresh-run solver.
using namespace Eigen;
using namespace std;

typedef std::complex<double> Complex;
const Complex I(0, 1);//defines i=sqrt(-1).

//defining different functions required for the simulation. DO NOT ALTER ANYTHING IN THE FOLLOWING PART UNTIL AND UNLESS YOU WANT TO CHANGE SOMETHING IN THE CODE FUNDAMENTALLY

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
complex<double> Velocity_calculation(const Ref<const VectorXcd> b, const Ref<const VectorXcd> b_new, int L, complex<double> V_i, double m, double g)
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

