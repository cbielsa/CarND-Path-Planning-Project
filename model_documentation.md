# Model documentation
## Source code files in this repository
* `spline.h` is a cubic spline interpolation library used by the Path Planner to fit trajectories `y(x)` in local car coordinates to a number of sample points.
The path planner only uses the constructor, member function `set_points` and `operator()`. Reference: http://kluge.in-chemnitz.de/opensource/spline/
* `main.cpp` contains the main program, as well the following auxiliary functions:
  * `ClosestWaypoing()` and `NextWaypoint()`, to calculate closest and next waypoints.
  * `getFrenet()` and `getXY()`, to transform between global cartesian and Frenet coordinates.
  Note that in addition to transforming the input position, these funtions also return a rotation matrix that allows to transform
  other vectors s.a. velocity and acceleration.
  * `lane_index()`, to calculate the lane index (`0`, `1` or `2`) from Frenet coordinate `d`.

## Algorithm outline
The high level steps of the algorithm are in `int main()` in `main.cpp`. At each time cycle:
* The telemetry from the simulator is extracted, containing the ego state (`x`, `y`, `s`, `d`, `yaw`, `speed`),
the ego previous path, the ego end path `s` and `d` Frenet coordinates, and the data from sensor fusion `sensor_fusion`
with the estimated position of each vehicle on the road within sensors range.
* The ego reference state (position, velocity and acceleration) is estimated.
If the previous path contains fewer than 2 points, then the reference state is set to the current ego state;
else, the reference state is set to the state at the end of the previous path.
Note that, if the previous path contains more than 2 points, then also the ego acceleration is estimated with divided differences.
This step is performed in section `[1] Estimate reference state` of `main()`.
* Target lane and speed are selected to maximize ego speed going forward, taking into account the position and speed of the other cars
on the road, extracted from `sensor_fusion`. A change of lane is only selected if its is both safe and advantageous to ego.
The target velocity is set to the maximum allowed velocity whenever the target lane has free space ahead, else it is set relative to
the velocity of the car ahead.
After this step, variables `target_d` and `target_vel` are set to the optimum values.
This step is performed in section `[2] Choose target lane and speed` of `main()`.
* The future trajectory is defined in local car xy-coordinates (origin at ego position, x-axis along Frenet `s`,
y-axis towards car left-side). Then a spline object `spline_traj` is fit to the trajectory points.
Note that `spline_traj` contains the trajectory, not the law of motion, i.e. it interpolates `y(x)` in local car coordinates but
does not contain information about ego position along the trajectory at a given time `t`. The law of motion is defined separetly in
the next step. To achieve trajectory continuity, the first two points in the trajectory correspond to the last two points of the
previously commanded path. Then, three more points are added 30, 60 and 90 meters ahead of the reference position, with Frenet
coordinate `target_d` (which was selected in the previous step).
This step is performed in section `[3] Define trajectory` of `main()`.
* Finally, the law of motion is defined, based on: velocity, acceleration and jerk at the reference state;
the trajectory `y(x)` computed in the previous step; and the target velocity `target_vel` selected before to maximize ego
speed going forward while keeping the ego vehicle safe. Note that the target velocity is reached as soon as physically possible
without exceeding maximum velocity, acceleration and jerk values.
This step is performed in section `[4] Define law of motion (trajectory as a function of time)` of `main()`.
In particular:
  * Two vectors `next_x_vals` and `next_y_vals` are initialized to the path points left from the previously commanded path
  (these vectors contain path point global xy-coordinates at `DT = 0.02 second` time steps).
  * New path points are appended to `next_x_vals` and `next_y_vals` at `DT` steps up to a size of `N = 50` is reached. For each
  new path point:
    * the local car x-coordinate (along Frenet `s`) is chosen for ego to reach the target velocity as quick as possible without
    exceeding maximum acceleration and jerk constraints.
    * the local car-y coordinate is obtained by spline interpolation in the previously generated trajectory: `y = spline_traj(x)`.
    * the path point is converted from car local to global xy coordinates and appended to `next_x_vals` and `next_y_vals`.
  * Once the path is complete, it is commanded to the simulator and the processing cycle is complete.
