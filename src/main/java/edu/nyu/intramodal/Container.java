package edu.nyu.intramodal;

import java.util.List;

public class Container {

    Integer pickup_index;
    Integer v1_id, v2_id;
    Integer v1_latest_arr_t, v2_submission_t;
    Integer number_of_transfer_nodes;
    Integer v1_transfer_waiting_t, v2_transfer_waiting_t;
    Double v1_cost;
    List<Vehicle.Route_stop> v1_route;
    Double v2_cost;
    List<Vehicle.Route_stop> v2_route;

    //used for without transfer result
    Container (Integer v_id, Double v_cost, List<Vehicle.Route_stop> v_route){
        this.v1_id = v_id;
        this.v1_cost = v_cost;
        this.v1_route = v_route;
    }

    //used for first vehicle of transfer
//    Container (Integer v_id, Integer pickup_index, Float v_cost, List<Vehicle.Route_stop> v_route){
//        this.v_id = v_id;
//        this.pickup_index = pickup_index;
//        this.v1_cost = v_cost;
//        this.v1_route = v_route;
//    }

    //used for transfer result
    Container (Integer v1_id, Integer v2_id, Double v1_cost, Double v2_cost, List<Vehicle.Route_stop> v1_route, List<Vehicle.Route_stop> v2_route,
               Integer v1_latest_arr_t, Integer v2_submission_t, Integer number_of_transfer_nodes, Integer v1_transfer_waiting_t,
               Integer v2_transfer_waiting_t){
        this.v1_id = v1_id;
        this.v2_id = v2_id;
        this.v1_cost = v1_cost;
        this.v2_cost = v2_cost;
        this.v1_route = v1_route;
        this.v2_route = v2_route;
        this.v1_latest_arr_t = v1_latest_arr_t;
        this.v2_submission_t = v2_submission_t;
        this.number_of_transfer_nodes = number_of_transfer_nodes;
        this.v1_transfer_waiting_t = v1_transfer_waiting_t;
        this.v2_transfer_waiting_t = v2_transfer_waiting_t;
    }
}
