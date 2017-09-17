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

## Algorithm outline
The high level steps of the algorithm are in `int main()` in `main.cpp`. At each time step:
* The telemetry from the simulator is extracted, containing the ego state (`x`, `y`, `s`, `d`, `yaw`, `speed`), the ego previous path, the ego end path `s` and `d` Frenet coordinates, and the data from sensor fusion `sensor_fusion` with the estimated position of each vehicle on the road within sensors range.
* The ego reference state (position, velocity and acceleration) is estimated. If the previous path contains fewer than 2 points, then the reference state is set to the current ego state; else, the reference state is set to the state at the end of the previous path. Note that, if the previous path contains more than 2 points, then also the ego acceleration is estimated with divided differences.
* Select the target state. This is the most complex step of the algorithm, consisting of the following steps:
 * Generate a list of candidate target states (after `N*DT` seconds):
  * `N_s` values of target s are considered, from `ref_s` to `ref_s + N*DT*max_speed`
  * `N_vs` values of target velocitiy along s are considered, from `0` to `max_speed`
  * `N_d` values of target d are considered (2., 6. and 10.) if lane change is enabled
 * For each candidate target state, generate minimum-jerk trajectory in Frenet coordinates (two order 5 polynomials `s(t)`, `d(t)`) from the reference and target states
 * Assign a cost to the trajectory to each candidate state. The cost function accounts for factors such as: close distance or impact to other vehicles in the road; excessive velocity, acceleration or jerk; and average speed and distance travelled
 * The target with the lowest cost is selected as the target state.
