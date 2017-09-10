
#include <vector>

using namespace std;


class Jmt
{

 public:

    //! Default constructor
    Jmt();

    //! Constructor from initial state, final state and time between states
    Jmt(
      const vector<double> &vState0,
      const vector<double> &vState1,
      double T );

    //! Calculate state at given input elapsed time
    double operator()( double dt ) const;

    //! Calculate first and second time derivatives of state at given elapsed time
    void velAccJerk( double dt, double &vel, double &acc, double &jerk ) const;

    //! Calculate max abs value of 1st and 2nd time derivatives for given
    //! input elapsed times
    void maxAbsVelAccJerk(
      const vector<double> &dt, double &maxVel, double &maxAcc, double &maxJerk ) const;

 private:

    //! Vector of coeffs of polynomial 
    //! s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5
    vector<double> _coeff;

};
