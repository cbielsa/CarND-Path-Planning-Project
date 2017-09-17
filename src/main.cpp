#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#include "jmt.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

double norm(double x, double y)
{
	return sqrt(x*x + y*y);
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/2)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}


// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
//   - the output vector contains the {s,d} of input position vector {x,y}
//   - xy2sd contains the transformation matrix from {x,y} to {s,d} coordinates
//     as a vector of size 4 (useful to rotate direction vectors)
vector<double> getFrenet(
	double x, double y, double theta,
	vector<double> maps_x, vector<double> maps_y,
	vector<double> &xy2sd )
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);
	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
		prev_wp  = maps_x.size()-1;

	double us_x = maps_x[next_wp]-maps_x[prev_wp];
	double us_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double us_norm2 = us_x*us_x+us_y*us_y;
	double proj_norm = (x_x*us_x+x_y*us_y)/us_norm2;
	double proj_x = proj_norm*us_x;
	double proj_y = proj_norm*us_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
		frenet_d *= -1;

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);

	frenet_s += distance(0,0,proj_x,proj_y);

	// normalize vector along s coordinate
	double us_norm = sqrt(us_norm2);
	us_x /= us_norm;
	us_y /= us_norm;

	// calculate rotation matrix to transform coords from XY to SD frame
	xy2sd = { us_x, us_y, us_y, -us_x };

	return {frenet_s, frenet_d};
}


// Transform from Frenet s,d coordinates to Cartesian x,y:
//   - the output vector contains the {x,y} of input position vector {s,d}
//   - sd2xy contains the transformation matrix from {s,d} to {x,d} coordinates
//     as a vector of size 4 (useful to rotate direction vectors)
vector<double> getXY(
	double s, double d,
	vector<double> maps_s, vector<double> maps_x, vector<double> maps_y,
	vector<double> &sd2xy )
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	// rotation matrix from sd to xy coords
	double cosh = cos(heading);
	double sinh = sin(heading);
	sd2xy = { cosh, -sinh, sinh, cosh };

	return {x, y};
}


//! Determine lane index for a given d value
int lane_index( double d )
{
	int i;

    if( d>0. && d<4. )
      	i = 0;
    else if( d < 8. )
      	i = 1;
    else if( d < 12. )
        i = 2;
    else
      	cout << "ERROR in lane_index: invalid d : " << d << endl;
      		
    return i;
}



// Main program ----------

int main()
{

  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }


  h.onMessage(
  	[ &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,
  	  &map_waypoints_dy ]
  	( uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode)
  	{

    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;

    // time-step [s]
  	const double DT = 0.02;
  	const double DT2 = DT*DT;

    // size of planned horizon
    int N = 50;  // N*DT set to 1 second

    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2')
    {

      auto s = hasData(data);

      if (s != "")
      {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry")
        {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = deg2rad(j[1]["yaw"]);     // [rad]
          	double car_speed = static_cast<double>(j[1]["speed"])*0.44704;  // [m/s]

          	// Previous path data given to the Planner
          	vector<double> previous_path_x = j[1]["previous_path_x"];
          	vector<double> previous_path_y = j[1]["previous_path_y"];
          	
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	

          	// [0] CONSTRAINT DEFINITION ==========================================

          	// max speed allowed
          	const double max_vel = 22.30-1.0;  // 22.30 m/s = 50 mph

          	// max acceleration allowed
          	const double max_acc = 10.-3.;  // m/s^2

          	// max jerk allowed
          	const double max_jerk = 90.;  // m/s^3



          	// [1] ESTIMATE REFERENCE STATE ======================================

			// set reference car state from which to resume path planning
          	// either at current car state or at last point of previous path

          	size_t prev_size = previous_path_x.size();

          	double ref_x, ref_y, ref_s, ref_d, ref_yaw;
          	double ref_vx, ref_vy, ref_vs, ref_vd;
          	double ref_x_prev, ref_y_prev;

          	// acceleration
          	double ref_ax = 0.;
          	double ref_ay = 0.;
          	double ref_as = 0.;
          	double ref_ad = 0.;


          	// if less than two points in previous path, set ref state to current state
          	if( prev_size < 2 )
          	{

          		// position
          		ref_x = car_x;
          		ref_y = car_y;
          		ref_s = car_s;
          		ref_d = car_d;
          		ref_yaw = car_yaw;  // [rad]

          		// velocity
         		ref_vx = car_speed*cos(car_yaw);
          		ref_vy = car_speed*sin(car_yaw);

          		// transformation matrix from xy to sd frame
				vector<double> xy2sd(4);
          		getFrenet(ref_x, ref_y, ref_yaw, map_waypoints_x, map_waypoints_y, xy2sd);
          		
          		// rotate velocity to Frenet frame
          		ref_vs = ref_vx*xy2sd[0] + ref_vy*xy2sd[1];
          		ref_vd = ref_vx*xy2sd[2] + ref_vy*xy2sd[3];

          		if( ref_vx > 0 )
          		{
          			ref_x_prev = ref_x-ref_vx*DT;
          			ref_y_prev = ref_y-ref_vy*DT;
          		}
          		else
          		{
          			ref_x_prev = ref_x-1.*DT*cos(car_yaw);
          			ref_y_prev = ref_y-1.*DT*sin(car_yaw);          			
          		}

          	}

          	// if at least two points, set ref state to last point of previously calc path
          	else
          	{

          		// position
          		ref_x = previous_path_x[prev_size-1];
          		ref_y = previous_path_y[prev_size-1];
          		ref_s = end_path_s;
          		ref_d = end_path_d;

          		ref_x_prev = previous_path_x[prev_size-2];
          		ref_y_prev = previous_path_y[prev_size-2];

          		ref_yaw = atan2( ref_y-ref_y_prev, ref_x-ref_x_prev );

          		// velocity
          		ref_vx = (ref_x - previous_path_x[prev_size-2])/DT;
          		ref_vy = (ref_y - previous_path_y[prev_size-2])/DT;

          		// if at least three points, calculate reference acceleration
          		if( prev_size > 2 )
          		{
          			ref_ax = (ref_x - 2*previous_path_x[prev_size-2]
          				            + previous_path_x[prev_size-3] )/DT2;
          			ref_ay = (ref_y - 2*previous_path_y[prev_size-2]
          				            + previous_path_y[prev_size-3] )/DT2;
          		}

          		// velocity and acceleration in Frenet coords
          		vector<double> xy2sd(4);
          		auto out = getFrenet(ref_x, ref_y, ref_yaw, map_waypoints_x, map_waypoints_y, xy2sd);

          	    ref_vs = xy2sd[0]*ref_vx + xy2sd[1]*ref_vy;
          	    ref_vd = xy2sd[2]*ref_vx + xy2sd[3]*ref_vy;
          	    ref_as = xy2sd[0]*ref_ax + xy2sd[1]*ref_ay;
          	    ref_ad = xy2sd[2]*ref_ax + xy2sd[3]*ref_ay;

          	}

          	if( ref_vs < 0. )
          	{
          		cout << endl;
          		cout << "negative ref_vs!" << endl;
          		cout << "*** END OF TRACK REACHED ***" << endl;
          		cout << endl;
          	}

       

          	// [2] CHOOSE TARGET LANE AND SPEED ===================================

		// sensor_fusion[i] : {id, x, y, vx, vy, s, d} with vx, vy in [m/s]

          	// s and speed of next car on each lane
          	vector<bool> next_car = {false, false, false};
          	vector<double> next_car_dist = {1e12, 1e12, 1e12};
          	vector<double> next_car_v = {0., 0., 0.};

          	vector<bool> prev_car = {false, false, false};
          	vector<double> prev_car_dist = {1e12, 1e12, 1e12};

    		// go over cars detected by sensor fusion system
    		for( size_t i=0; i<sensor_fusion.size(); ++i )
    		{

	        	// obstacle car d (we assume it remains in same lane)
	        	double obs_car_d = sensor_fusion[i][6];

	        	// propagate obstacle car s until end of previous path time,
        		// assuming constant speed
        		double obs_car_vx = sensor_fusion[i][3];
        		double obs_car_vy = sensor_fusion[i][4];
        		double obs_car_speed = norm(obs_car_vx, obs_car_vy);
      			double obs_car_s = sensor_fusion[i][5];
      			obs_car_s += static_cast<double>(prev_size)*DT*obs_car_speed;

      			// distance to ego at reference time (end of prev path)
      			double obs_car_ds = obs_car_s - ref_s;		
      			int obs_car_lane = lane_index(obs_car_d);

      			// case obstacle car in front
      			if( obs_car_ds > 0. )
      			{
      				next_car[obs_car_lane] = true;

      				if( obs_car_ds < next_car_dist[obs_car_lane] )
      				{
      					next_car_dist[obs_car_lane] = obs_car_ds;
      					next_car_v[obs_car_lane] = obs_car_speed;
      				}
      			}

      			// case obstacle car behind
      			else
      			{

      				prev_car[obs_car_lane] = true;

     				if( -obs_car_ds < prev_car_dist[obs_car_lane] )
      				{
      					prev_car_dist[obs_car_lane] = -obs_car_ds;
      				}
      			}

      		}

		
      		// initialise target lane to lane in reference state
      		// and target velocity to maximum velocity
      		int target_lane = lane_index(ref_d);
      		double target_vel = max_vel;
      		bool bTargetSelected = false;

      		// case no car in front: keep lane and target max velocity
      		if( !next_car[target_lane] || next_car_dist[target_lane] > 50. )
      		{
      			target_lane = target_lane;
      			target_vel = max_vel;
      			bTargetSelected = true;

      			cout << "no car ahead" << endl;
      		}

      		// case car in front while in left or right lane: consider middle lane
      		else if( target_lane == 0 || target_lane == 2 )
      		{
      			// change lane if it is safe and car in candidate lane is faster than car in current lane
      			int candidate_lane = 1;

      			if( !next_car[candidate_lane] ||
      				( next_car_dist[candidate_lane] > next_car_dist[target_lane]
      			        && next_car_v[candidate_lane] > next_car_v[target_lane] ) 
      				&& ( !prev_car[candidate_lane] || prev_car_dist[candidate_lane]>35. ) )
      			{
      				target_lane = candidate_lane;
      				if( next_car[candidate_lane] )
      					target_vel = next_car_v[target_lane] - 0.5;
      				else
      					target_vel = max_vel;
      				bTargetSelected = true;

      				cout << "target middle lane" << endl;
      			}
      			    		
      		}

      		// case car in front while in middle lane: consider left and right lanes
      		else if( target_lane == 1 )
      		{
      			// consider right if it is safe and car in candidate lane is faster than car in current lane
      			int candidate_lane = 2;

      			if( !next_car[candidate_lane] ||
      				( next_car_dist[candidate_lane] > next_car_dist[target_lane]
      			        && next_car_v[candidate_lane] > next_car_v[target_lane] ) 
      				&& ( !prev_car[candidate_lane] || prev_car_dist[candidate_lane]>35. ) )
      			{
      				target_lane = candidate_lane;
      				if( next_car[candidate_lane] )
      					target_vel = next_car_v[target_lane] - 0.5;
      				else
      					target_vel = max_vel;
      				bTargetSelected = true;

      				cout << "target right lane" << endl;
      			}

      			// consider left if it is safe and car in candidate lane is faster than car in current lane
      			candidate_lane = 0;

      			if( !next_car[candidate_lane] ||
      				( next_car_dist[candidate_lane] > next_car_dist[target_lane]
      			        && next_car_v[candidate_lane] > next_car_v[target_lane] ) 
      				&& ( !prev_car[candidate_lane] || prev_car_dist[candidate_lane]>35. ) )
      			{
      				target_lane = candidate_lane;
      				if( next_car[candidate_lane] )
      					target_vel = next_car_v[target_lane] - 0.5;
      				else
      					target_vel = max_vel;
      				bTargetSelected = true;

      				cout << "target left lane" << endl;
      			}
      		}

      		// case car in front but lane change is not safe or advantageous
      		// stay in middle lane but reduce speed
      		if( !bTargetSelected )
      		{
      			target_vel = next_car_v[target_lane] - 0.5;

      			cout << "keep lane and reduce speed" << endl;
      		}

      		//cout << "target_lane : " << target_lane << endl;

      		double target_d = 2. + 4.*target_lane;



      		// [3] DEFINE TRAJECTORY =================================================

          	// trajectory points in global xy coordinates
          	vector<double> ptsx, ptsy;

          	// initialize next path with two points: at time step before reference
          	// and at reference time, to ensure continuity with previously commanded path
          	ptsx.push_back(ref_x_prev);
          	ptsx.push_back(ref_x);

          	ptsy.push_back(ref_y_prev);
          	ptsy.push_back(ref_y);

          	// add three target waypoints 30, 60 and 90 m ahead of reference state
			vector<double> sd2xy(4);
			auto xy = getXY(
          		ref_s+30., target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y, sd2xy);
			ptsx.push_back(xy[0]); ptsy.push_back(xy[1]);

			xy = getXY(
          		ref_s+60., target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y, sd2xy);
			ptsx.push_back(xy[0]); ptsy.push_back(xy[1]);

			xy = getXY(
          		ref_s+90., target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y, sd2xy);
			ptsx.push_back(xy[0]); ptsy.push_back(xy[1]);

			// transform trajectory points from global to local car xy frame (x along car velocity, y towards car left)
			for( size_t i=0; i<ptsx.size(); ++i )
			{

				double shift_x = ptsx[i]-ref_x;
				double shift_y = ptsy[i]-ref_y;

				ptsx[i] =  shift_x*cos(ref_yaw) + shift_y*sin(ref_yaw);
				ptsy[i] = -shift_x*sin(ref_yaw) + shift_y*cos(ref_yaw);
			}

			// create a spline and fit it to the trajectory, in local car xy frame
			tk::spline spline_traj;
			spline_traj.set_points(ptsx, ptsy);



          	// [4] DEFINE LAW OF MOTION (TRAJECTORY AS FUNCTION OF TIME) ==========

         	// Sample trajectory at DT steps

     		// copy waypoints of previously calculated path to this path
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	for( size_t i=0; i<prev_size; ++i )
          	{
          		next_x_vals.push_back( previous_path_x[i] );
          		next_y_vals.push_back( previous_path_y[i] );
          	}

          	// start position, velocity and acceleration
          	double s = 0.;
          	double v = ref_vs;
          	double a = ref_as;

          	// brake or accelerate ?
          	bool bAccel = ( target_vel > v ? true : false );

          	// start jerk
          	double j = ( bAccel ? 1 : -1  )*max_jerk;


          	//cout << "-------------------" << endl;
          	//cout << "accel? " << bAccel << ", target_vel : " << target_vel << ", ref_vs : " << ref_vs << ", ref_as : " << ref_as << endl;

          	// define quickest motion to reach target velocity w/o kinematic violations
          	for( size_t i=0; i<N-prev_size; ++i )
          	{

          		// next position
          		s += v*DT;

          		// next velocity
          		v += a*DT;
          		
          		// case target velocity reached
          		if(    (bAccel && (v > target_vel))
          			|| ((!bAccel) && (v < target_vel)) )
          		{
          			v = target_vel;
          		    a = 0.;
          		    j = 0.;
          		}

          		// next acceleration
          		a += j*DT;
          		if( a > max_acc )
          		{
          			a = max_acc;
          			if( bAccel )
          				j = 0.;
          		}
          		else if( a < -max_acc )
          		{
          			a = -max_acc;
          			if( !bAccel )
          				j = 0.;
          		}

          		//cout << "j:" << j << " a:" << a << " v:" << v << " s:" << s <<endl;

          		// assume road stretch from reference position to end of path is straight
          		double x = s;
          		double y = spline_traj(x);

          		// transform path point from local to global coordinates
          		double xg = ref_x + x*cos(ref_yaw) - y*sin(ref_yaw);
          		double yg = ref_y + x*sin(ref_yaw) + y*cos(ref_yaw);

          		// add point to path
          		next_x_vals.push_back(xg);
          		next_y_vals.push_back(yg);

          	}


          	// pass to simulator path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT); 
        }
      
      }

      else
      {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    
    }
  });


  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });


  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });


  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });


  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
