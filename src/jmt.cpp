
#include "jmt.h"

//#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using Eigen::VectorXd;
using Eigen::MatrixXd;


Jmt::Jmt()
  : _coeff(6)
{}


Jmt::Jmt( const vector<double> &vState0, const vector<double> &vState1, double T )
  : _coeff(6)
{

  // parse input arguments
  const double& s0 = vState0[0];
  const double& v0 = vState0[1];
  const double& a0 = vState0[2];
  const double& s1 = vState1[0];
  const double& v1 = vState1[1];
  const double& a1 = vState1[2];
  
  // auxiliary variables
  double T2 = T*T;
  double T3 = T2*T;
  double T4 = T2*T2;
  double T5 = T2*T3;
      
  // define linear system for last three coefficients
  VectorXd b(3);
  b << s1-s0-v0*T-a0*T2/2.,
       v1-v0-a0*T,
       a1-a0;
  
  MatrixXd A(3,3);
  A <<    T3,     T4,     T5,
       3.*T2,  4.*T3,  5.*T4,
       6.*T,  12.*T2, 20.*T3;
  
  // solve linear system for last three coefficients
  VectorXd x = A.colPivHouseholderQr().solve(b);

  // set elements of coefficient vector
  _coeff[0] = s0;
  _coeff[1] = v0;
  _coeff[2] = a0/2.;
  _coeff[3] = x[0];
  _coeff[4] = x[1];
  _coeff[5] = x[2];    
}


double Jmt::operator()( double dt ) const
{
  double dt2 = dt*dt;
  double dt3 = dt2*dt;
  double dt4 = dt2*dt2;
  return  _coeff[0] + _coeff[1]*dt + _coeff[2]*dt2 + _coeff[3]*dt3
                    + _coeff[4]*dt4 + _coeff[5]*dt4*dt;
}


void Jmt::velAccJerk( double dt, double &vel, double &acc, double &jerk ) const
{
  double dt2 = dt*dt;
  double dt3 = dt2*dt;
  double dt4 = dt2*dt2;
  
  vel = _coeff[1] + 2*_coeff[2]*dt + 3*_coeff[3]*dt2
        + 4*_coeff[4]*dt3 + 5*_coeff[5]*dt4;
  acc = 2*_coeff[2] + 6*_coeff[3]*dt + 12*_coeff[4]*dt2
        + 20*_coeff[5]*dt3;
  jerk = 6*_coeff[3] + 24*_coeff[4]*dt + 60*_coeff[5]*dt2;
  return;
}


void Jmt::maxAbsVelAccJerk(
  const vector<double> &dt, double &maxVel, double &maxAcc, double &maxJerk ) const
{
  maxVel = 0.;
  maxAcc = 0.;
  maxJerk = 0.;

  for( size_t i=0; i<dt.size(); ++i )
  {
    double vel, acc, jerk;
    velAccJerk(dt[i], vel, acc, jerk);
    if( fabs(vel) > maxVel )
      maxVel = vel;
    if( fabs(acc) > maxAcc )
      maxAcc = acc;
    if( fabs(jerk) > maxJerk )
      maxJerk = jerk;
  }

  return;
}
