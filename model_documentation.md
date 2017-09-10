# Model documentation
## Source code files in this repository
* `jmt.h` and `jmt.cpp` declare and define a class to calculate a minimum-jerk polynomial from start and end states and the time duration
between those states, where a state is given by position, velocity and acceleration along a given coordinate. After construction of the
polynomial, the class provides member functions to calculate position (`operator()()`), velocity, acceleration and jerk (`velAccJerk()`)
at a given input time, as well as the maximum velocity, acceleration and jerk amongst those calculated at given input times (`maxAbsVelAccJerk()`).
* `main.cpp` contains the main program, as well the following auxiliary functions:
  * `ClosestWaypoing()` and `NextWaypoint()`, to calculate closest and next waypoints.
  * `getFrenet()` and `getXY()`, to transform between global cartesian and Frenet coordinates.
  Note that in addition to transforming the input position, these funtions also return a rotation matrix that allows to transform
  other vectors s.a. velocity and acceleration.
  * `costFunction()`, to calculate the cost function for a trajectory with given initial and final states. For the given states,
  a minimum-jerk trajectory is constructed in Frenet coordinates, and the following factors contribute to the cost:
    * excessive speed along `s` at any point along the candidate trajectory
    * excessive acceleration (both `s` and `d` coordinates are considered)
    * excessive jerk (both `s` and `d` coordinates are considered)
    * low target speend along `s`
    * low target `s` (which results in low average speed)
    * changes of lane (to prevent unnecessary lane changes)
    * deviations from central lane (to prefer central lane `1`, since it gives more options to select lane change in the future)
    * obstacle cars with similar `d` than ego and a near `s` (risk of collision)
