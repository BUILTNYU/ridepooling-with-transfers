package edu.nyu.intramodal;

import java.util.ArrayList;
import java.util.List;

public class Vehicle {

    int current_loc, current_load, tt_to_next_node;
    double cost;
    List<Route_stop> route_stops;
    List<Integer> active_requests;
    List<Integer> route_nodes;

    Vehicle(int current_loc, int current_load, int tt_to_next_node, double cost, List<Route_stop> route_stops,
            List<Integer> route_nodes, List<Integer> active_requests){

        this.current_loc = current_loc;
        this.current_load = current_load;
        this.tt_to_next_node = tt_to_next_node;
        this.cost = cost;
        this.route_stops = route_stops;
        this.route_nodes = route_nodes;
        this.active_requests = active_requests;
    }

    static class Route_stop implements Cloneable{
        Integer stop_id, exp_stop_arr_t;
        List<Integer> pickup_ids, dropoff_ids;

        Route_stop(Integer stop_id, Integer exp_stop_arr_t, List<Integer> pickup_ids, List<Integer> dropoff_ids){
            this.stop_id = stop_id;
            this.exp_stop_arr_t = exp_stop_arr_t;
            this.pickup_ids = pickup_ids;
            this.dropoff_ids = dropoff_ids;
        }

        public Route_stop clone() throws CloneNotSupportedException {
            Route_stop copy = (Route_stop) super.clone();
            copy.pickup_ids = new ArrayList<>(pickup_ids);
            copy.dropoff_ids = new ArrayList<>(dropoff_ids);
            return copy;
        }
    }
}
