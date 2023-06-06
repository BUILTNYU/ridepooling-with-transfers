package edu.nyu.intramodal;

import java.util.HashMap;
import java.util.List;

public class Data {

    HashMap<Integer, VehicleSim> vehicles;

    Data(HashMap<Integer, VehicleSim> vehicles_sim){
        this.vehicles = vehicles_sim;
    }

    static class VehicleSim {
        Integer n_visited_stops, n_served_pax, total_idle_t, total_distance, initial_loc, total_empty_distance;
        List<Vehicle.Route_stop> route_stops;

        VehicleSim (Integer n_visited_stops, Integer n_served_pax, Integer total_idle_t, List<Vehicle.Route_stop> route_stops, Integer total_distance, Integer initial_loc, Integer total_empty_distance){
            //TODO remove n_visited_stops and n_served_pax?
            this.n_visited_stops = n_visited_stops;
            this.n_served_pax = n_served_pax;
            this.total_idle_t = total_idle_t;
            this.route_stops = route_stops;
            this.total_distance = total_distance;
            this.initial_loc = initial_loc;
            this.total_empty_distance = total_empty_distance;
        }
    }

//    static class RequestSim {
//        Integer pickup_t, dropoff_t;
//        RequestSim (Integer pickup_t, Integer dropoff_t){
//            this.pickup_t = pickup_t;
//            this.dropoff_t = dropoff_t;
//        }
//    }
}
