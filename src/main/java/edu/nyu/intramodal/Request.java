package edu.nyu.intramodal;

import java.util.ArrayList;
import java.util.List;

public class Request {

    int origin, dest, submission_t, p_latest_arr_t;
    Integer pickup_t, dropoff_t, assigned_veh;
    boolean transfer;
    String status;

    Request(int origin, int dest, int submission_t, int p_latest_arr_t, Integer pickup_t, Integer dropoff_t,
            String status, List<Integer> assigned_vehs, boolean transfer){
        this.origin = origin;
        this.dest = dest;
        this. submission_t = submission_t;
        this.p_latest_arr_t = p_latest_arr_t;
        this.pickup_t = pickup_t;
        this.dropoff_t = dropoff_t;
        this.status = status;
        this.assigned_veh = assigned_veh;
        this.transfer = transfer;
    }
}
