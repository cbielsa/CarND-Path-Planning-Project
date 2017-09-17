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

	//cout << "heading error to next wy [deg] : " << rad2deg(angle) << endl; 

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

	//cout << "next waypoint : " << next_wp << endl;

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



// cost function to evaluate a given candidate trajectory

double costFunction(
	// reference (initial) state
    double ref_s, double ref_vs, double ref_as,
    double ref_d, double ref_vd, double ref_ad,
    // target (final) state
    double s, double vs, double as,
    double d, double vd, double ad,
    // thresholds
    double max_speed, double max_acc2, double max_jerk2,
    // input from sensor fusion
    const vector<vector<double> > &sensor_fusion,
    // vector of delta-times for threshold check
    const vector<double> &vDeltaTimes,
    // number of prediction horizon steps and delta-time [s]
    size_t N, double DT,
    // size of previous prediction horizon
    size_t prev_size,
    // debug info
    bool bDebugInfo, string &sDebugMsg )
{

    // initialise cost
    double cost = 0.;
    sDebugMsg = "";

    // calculate minimum-jerk profile for s coord
    Jmt jmtS(
    	{ref_s, ref_vs, ref_as},
        {s, vs, as}, N*DT );

    // calculate max abs of vel and acc along s
    double maxVel, maxAccS, maxJerkS;
    jmtS.maxAbsVelAccJerk(
  		vDeltaTimes, maxVel, maxAccS, maxJerkS);

	// penalize excessive speed along s
	if( maxVel > max_speed )
	{
		cost += 1000;
		if( bDebugInfo )
			sDebugMsg += "speed exceeded, ";
	}

	// penalise low target speed along s
	cost += fabs( max_speed - vs )/max_speed;

	// penalise low average speed along s
	double max_s = ref_s + max_speed*N*DT;
	cost += fabs( max_s - s )/max_s;

	// calculate minimum-jerk profile for d coord
    Jmt jmtD(
        {ref_d, ref_vd, ref_ad},
        {d, vd, ad}, N*DT );

    // calculate max abs of vel and acc along d
    double maxAccD, maxJerkD;
	jmtD.maxAbsVelAccJerk(
  		vDeltaTimes, maxVel, maxAccD, maxJerkD );

	// penalize excessive speed along d
	//if( maxVel > max_speed )
	//	cost += 1000;

	// penalize excessive acceleration
	if( maxAccS*maxAccS + maxAccD*maxAccD > max_acc2 )
	{
		cost += 1000;
		if( bDebugInfo )
			sDebugMsg += "acc exceeded, ";
	}

	/*
	// penalize excessive jerk
	if( maxJerkS*maxJerkS + maxJerkD*maxJerkD > max_jerk2 )
	{
		cost += 1000;
		if( bDebugInfo )
		{
			sDebugMsg += "jerk exceeded, ";
			if( maxJerkS > 8. )
				sDebugMsg += "large jerkS, ";
			if( maxJerkD > 8. )
				sDebugMsg += "large jerkD, ";
		}
	}
	*/

	// penalise lane changes
	cost += 0.2 * fabs(ref_d - d)/8.;

	// prefer middle lane
	if( fabs(d-6.) > 1.5 )
		cost += 0.15;
 
	// penalize being out of lane
	//cost += 10* fmod(d-2.,4) /2.;


	// penalize collisions -------------

	// sensor_fusion[i] : {id, x, y, vx, vy, s, d} with vx, vy in [m/s]

    // go over cars detected by sensor fusion system
    for( size_t i=0; i<sensor_fusion.size(); ++i )
    {

        // obstacle car d (we assume it remains in same lane)
        double car_d = sensor_fusion[i][6];

        // propagate obstacle car s until end of previous path time,
        // assuming consant speed
        double car_vx = sensor_fusion[i][3];
        double car_vy = sensor_fusion[i][4];
        double car_speed = sqrt(car_vx*car_vx+car_vy*car_vy);

      	double car_s = sensor_fusion[i][5];		
      	car_s += static_cast<double>(prev_size)*DT*car_speed;

       	// calculate distance to ego
       	double ds = car_s - s;
       	double dd = car_d - d;

       	//std::cout << "d : " << d << ", car_d : " << car_d << std::endl;

       	// if close to ego along d

       	if( fabs(dd) < 1. )
       	{
       		// penalize coming close to car in front or behind
       		if( ds > 0. )
       		{
     			cost += 1e5/(ds*ds*ds);
				if( bDebugInfo && ds<30. )
					sDebugMsg += "risk of collision, ";
     		}

/*
       		if( (ds > 0. && ds < 40.) || ( ds < 0. && -ds < 5. ) )
       		{
     			cost += 2000.;
				if( bDebugInfo )
					sDebugMsg += "risk of collision, ";
     		}
*/
       	}
      	//else if( fabs(dd) < 3. )
       	//{
       	//	// penalize coming close to car in front or behind
       	//	if( (ds > 0. && ds < 30.) || ( ds < 0. && -ds < 10. ) )
       	//		cost += 10.;
       	//}

    }

	return cost;
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

  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

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
    
    // time-step [s]
  	const double DT = 0.02;
  	const double DT2 = DT*DT;

    // size of planned horizon
    int N = 75;  // N*DT set to 1.5 second

  	// vector of delta-times later used to calculate max expected vel and acc
  	vector<double> vDeltaTimes(N+1);
  	for( size_t i=0; i<=N; ++i )
  		vDeltaTimes[i] = i*DT;


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

          	

          	// GUIDANCE =============================================

          	// max speed allowed
          	// (1 m/s margin to cope with simulation latency)
          	const double max_speed = 22.30;  // 22.30 m/s = 50 mph

          	// max acceleration allowed
          	// (2 m/s^2 margin to cope with simulation latency)
          	const double max_acc2 = pow(10.-1., 2);  // m/s^2

          	// max jerk allowed
          	// (1 m/s^3 margin to cope with simulation latency)
          	const double max_jerk2 = pow(10.-1., 2);  // m/s^3



          	// Estimate reference state ======================================

			// set reference car state from which to resume path planning
          	// either at current car state or at last point of previous path

          	size_t prev_size = previous_path_x.size();

          	double ref_x, ref_y, ref_s, ref_d, ref_yaw;
          	double ref_vx, ref_vy, ref_vs, ref_vd;

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
          	}

          	// if at least two points, set ref state to last point of previously calc path
          	else
          	{
          		// position
          		ref_x = previous_path_x[prev_size-1];
          		ref_y = previous_path_y[prev_size-1];
          		ref_s = end_path_s;
          		ref_d = end_path_d;

          		ref_yaw = atan2( ref_y-previous_path_y[prev_size-2],
          			             ref_x-previous_path_x[prev_size-2] );

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
          		getFrenet(ref_x, ref_y, ref_yaw, map_waypoints_x, map_waypoints_y, xy2sd);
          	    ref_vs = xy2sd[0]*ref_vx + xy2sd[1]*ref_vy;
          	    ref_vd = xy2sd[2]*ref_vx + xy2sd[3]*ref_vy;
          	    ref_as = xy2sd[0]*ref_ax + xy2sd[1]*ref_ay;
          	    ref_ad = xy2sd[2]*ref_ax + xy2sd[3]*ref_ay;
          	}


          	// Set target state ================================================

          	// generate list of candidate target states

          	const size_t N_s  = 12;
          	const size_t N_vs = 12;
          	const size_t N_d = 1;

          	vector<double> v_candidate_target_s(N_s);
          	vector<double> v_candidate_target_vs(N_vs);
          	vector<double> v_candidate_target_d(N_d);

          	for( size_t i=0; i<N_s; ++i )
          	{
          		v_candidate_target_s[i]
          			= ref_s +  (1.+static_cast<double>(i))/N_s *N*DT*max_speed;
          	}

          	for( size_t i=0; i<N_vs; ++i )
          	{
          		v_candidate_target_vs[i]
          			= (1.+static_cast<double>(i))/N_vs *max_speed;
          	}

          	//v_candidate_target_d[0] = 2.;   // left lane
          	//v_candidate_target_d[1] = 4.;
          	v_candidate_target_d[0] = 6.;   // middle lane
          	//v_candidate_target_d[3] = 8.;  
          	//v_candidate_target_d[2] = 10.;  // right lane


          	// for each target candidate, calculate min-jerk trajectory and cost
          	// select target candidate with lowest cost

          	double target_s, target_vs;
          	double target_as = 0;

          	double target_d;
          	double target_vd = 0.;
          	double target_ad = 0.;

          	double minCost = 10e12;

          	string sDebugMsg, sMsg;


          	for( size_t iS=0; iS<N_s; ++iS )
          	{
          		for( size_t iVS=0; iVS<N_vs; ++iVS )
          		{
          			for( size_t iD=0; iD<N_d; ++iD )
          			{

          				double cost = costFunction(
							// reference (initial) state
    						ref_s, ref_vs, ref_as, ref_d, ref_vd, ref_ad,
    						// target (final) state
    						v_candidate_target_s[iS], v_candidate_target_vs[iVS], target_as,
    						v_candidate_target_d[iD], target_vd, target_ad,
    						// thresholds
    						max_speed, max_acc2, max_jerk2,
   							// input from sensor fusion
    						sensor_fusion,
    						// vector of delta-times for threshold check
    						vDeltaTimes,
    						// number of prediction horizon steps and delta-time [s]
							N, DT,
							// size of previous prediction horizon
    						prev_size,
							// debug info
    						true, sDebugMsg );

			    		// check if this is best candidate so far
			    		if( cost < minCost )
			    		{
			    			minCost = cost;
			    			target_s = v_candidate_target_s[iS];
			    			target_vs = v_candidate_target_vs[iVS];
			    			target_d = v_candidate_target_d[iD];
			    			sMsg = sDebugMsg;
			    		}

          			}
          		}
          	}

          	std::cout << "min cost : " << minCost << ", msg : " << sMsg << std::endl;


          	// convert target state to global XY coordinates
          	vector<double> sd2xy(4);
          	auto out = getXY(
          		target_s, target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y, sd2xy);
          	double target_x = out[0];
          	double target_y = out[1];

          	double target_vx = sd2xy[0]*target_vs + sd2xy[1]*target_vd;
          	double target_vy = sd2xy[2]*target_vs + sd2xy[3]*target_vd;

          	double target_ax = 0.;
          	double target_ay = 0.;


          	// Calculate Minimum Jerk Trajectory for reference and target states =======

          	Jmt jmtX({ref_x, ref_vx, ref_ax}, {target_x, target_vx, target_ax}, N*DT);
          	Jmt jmtY({ref_y, ref_vy, ref_ay}, {target_y, target_vy, target_ay}, N*DT);

          	// Sample trajectory at DT steps

     		// copy waypoints of previously calculated path to this path
          	vector<double> next_x_vals = previous_path_x;
          	vector<double> next_y_vals = previous_path_y;

          	// extend this path with samples from minimim jerk trajectory
          	for( size_t i=0; i<N-prev_size; ++i )
          	{
          		next_x_vals.push_back( jmtX( (i+1)*DT ) );
          		next_y_vals.push_back( jmtY( (i+1)*DT ) );
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
