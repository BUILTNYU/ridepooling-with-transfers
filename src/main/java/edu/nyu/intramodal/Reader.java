package edu.nyu.intramodal;

import org.osgeo.proj4j.BasicCoordinateTransform;
import org.osgeo.proj4j.CRSFactory;
import org.osgeo.proj4j.CoordinateReferenceSystem;
import org.osgeo.proj4j.ProjCoordinate;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Reader {

    public static void read_network(Network network, String filename) throws IOException {

        //csv data structure: origin - destination - travel time - shortest path
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        br.readLine();
        while((line = br.readLine()) != null){
            String[] nn_data = line.split(",");
//                String[] s_path = nn_data[3].split(";");
//                ArrayList<Integer> shortest_path = new ArrayList<>();
//                for (int j=0; j<s_path.length; j++){
//                    shortest_path.add(Integer.parseInt(s_path[j]));
//                }
                network.addShortestPath(Integer.parseInt(nn_data[0]), Integer.parseInt(nn_data[1]), Integer.parseInt(nn_data[2]));
//            System.out.println("shortest path from node : "+ Integer.parseInt(nn_data[0]) + " to node " + Integer.parseInt(nn_data[1]) + " is : " + network.shortest_path[Integer.parseInt(nn_data[0])][Integer.parseInt(nn_data[1])]);
//            System.out.println("travel time from node : "+ Integer.parseInt(nn_data[0]) + " to node " + Integer.parseInt(nn_data[1]) + " is : " + network.dist[Integer.parseInt(nn_data[0])][Integer.parseInt(nn_data[1])]);
        }
    }

    public static void read_nodes (HashMap<Integer, Node> nodes, String filename) throws IOException {
        //csv data structure: network node id - X - Y - mapped node id
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        br.readLine();
        while((line = br.readLine()) != null){
            String[] line_data = line.split(",");
                nodes.put(Integer.parseInt(line_data[3]), new Node(Double.parseDouble(line_data[1]), Double.parseDouble(line_data[2]), Long.parseLong(line_data[0])));
        }
    }

    public static void read_stops (ArrayList<Integer> stop_ids, HashMap<Integer, Node> nodes, HashMap<Integer, Node> stops, String filename) throws IOException {
        //csv data structure: stop id - X - Y - ...
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        br.readLine();
        while((line = br.readLine()) != null){

            HashMap<Integer, double[]> nearest_node;
            String[] line_data = line.split(",");

            if (line_data[0].contains("r") == false){

                int former_stop_id = Integer.parseInt(line_data[0]);

                nearest_node = find_nearest_node(nodes, Double.parseDouble(line_data[1]), Double.parseDouble(line_data[2]));

                for (Map.Entry<Integer, double[]> item: nearest_node.entrySet()) {
                    stops.put(item.getKey(), new Node(item.getValue()[0], item.getValue()[1], former_stop_id));
                    if (stop_ids.contains(Integer.valueOf(item.getKey()))){
                        continue;
                    }
                    else {
                        stop_ids.add(item.getKey());
                    }
                }
            }
        }
    }

    public static void read_transfers (ArrayList<Integer> transfer_ids, HashMap<Integer, Node> stops, String filename) throws IOException {
        //csv data structure: stop id
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        while((line = br.readLine()) != null) {

            String[] line_data = line.split(",");

            if (line_data[0].contains("r") == false) {

                int former_stop_id = Integer.parseInt(line_data[0]);
                for (Map.Entry<Integer, Node> stop : stops.entrySet()){
                    if (stop.getValue().former_id == former_stop_id){
                        if (!transfer_ids.contains(Integer.valueOf(stop.getKey()))) {
                            transfer_ids.add(stop.getKey());
                        }
                        continue;
                    }
                }
            }
        }
    }

    public static void read_requests (HashMap<Integer, Request> requests, HashMap<Integer, Node> stops, Network network, String filename,
                                      int max_wait_t, float max_tt_p, int max_tt_min, int max_tt_add) throws IOException {

        //csv data structure: origin lng - origin lat - dest lng - dest lat - submission time
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        br.readLine();
        int d = 0;

        while((line = br.readLine()) != null){

            String[] line_data = line.split(",");
            HashMap<Integer, double[]> origin_stop;
            HashMap<Integer, double[]> dest_stop;
            Integer origin_stop_id = null;
            Integer dest_stop_id = null;

            ProjCoordinate origin = convert_coordinates(Double.parseDouble(line_data[0]), Double.parseDouble(line_data[1]));
            ProjCoordinate dest = convert_coordinates(Double.parseDouble(line_data[2]), Double.parseDouble(line_data[3]));

            origin_stop = find_nearest_node(stops, origin.x, origin.y);
            dest_stop = find_nearest_node(stops, dest.x, dest.y);

            for (Map.Entry<Integer, double[]> item: origin_stop.entrySet()) {
                origin_stop_id = item.getKey();
            }
            for (Map.Entry<Integer, double[]> item: dest_stop.entrySet()) {
                dest_stop_id = item.getKey();
            }

            //submission time
            int sub_time;
            try {
                sub_time = (int) Double.parseDouble(line_data[4]);
                while (sub_time > 86400){
                    sub_time -= 86400;
                }
            }
            catch (NumberFormatException e) {
                String pattern = "yyyy-MM-dd'T'HH:mm:ss'Z'";
                DateTimeFormatter formatter = DateTimeFormatter.ofPattern(pattern);
                LocalDateTime local_date_time = LocalDateTime.parse(line_data[4], formatter);
                sub_time = local_date_time.getHour() * 3600 + local_date_time.getMinute() * 60 + local_date_time.getSecond();
            }

            //adding the request to the requests object
            requests.put(d, load_request(origin_stop_id, dest_stop_id, sub_time, network.dist, max_wait_t, max_tt_p,
                    max_tt_min, max_tt_add));
            d += 1;
        }
    }

    public static void read_hub_locations (ArrayList<Integer> hub_locations, HashMap<Integer, Node> stops, String filename) throws IOException {
        //csv data structure: hub lng - hub lat
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        br.readLine();

        while((line = br.readLine()) != null){

            String[] line_data = line.split(",");
            HashMap<Integer, double[]> hub_stop;

            ProjCoordinate hub = convert_coordinates(Double.parseDouble(line_data[0]), Double.parseDouble(line_data[1]));

            hub_stop = find_nearest_node(stops, hub.x, hub.y);

            for (Map.Entry<Integer, double[]> item: hub_stop.entrySet()) {
                if (!hub_locations.contains(item.getKey())) {
                    hub_locations.add(item.getKey());
                }
            }
        }
    }

    public static double find_distance (double x1, double y1, double x2, double y2){
        double distance = Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2));
        return distance;
    }

    public static ProjCoordinate convert_coordinates (double lng, double lat) {
        //converting coordinates to EPSG:4326
        CRSFactory factory = new CRSFactory();
        CoordinateReferenceSystem srcCrs = factory.createFromName("EPSG:4326"); //it's equivalent to WGS84
        CoordinateReferenceSystem dstCrs = factory.createFromName("EPSG:25832");

        BasicCoordinateTransform transform = new BasicCoordinateTransform(srcCrs, dstCrs);

        // Note these are x, y so lng, lat
        ProjCoordinate srcCoord = new ProjCoordinate(lng, lat);
        ProjCoordinate dstCoord = new ProjCoordinate();

        // Writes result into dstCoord
        transform.transform(srcCoord, dstCoord);

        return dstCoord;
    }

    public static HashMap<Integer, double[]> find_nearest_node (HashMap<Integer, Node> nodes, double x, double y) {
        double shortest_dist = Math.pow(10, 10);
        double distance;
        HashMap<Integer, double[]> nearest_node = new HashMap<>();
        for (Map.Entry<Integer, Node> node: nodes.entrySet()){
            distance = find_distance(x, y, node.getValue().x, node.getValue().y);
            if (distance < shortest_dist){
                nearest_node.clear();
                nearest_node.put(Integer.valueOf(node.getKey()), new double[] {Double.valueOf(node.getValue().x), Double.valueOf(node.getValue().y)});
                shortest_dist = distance;
            }
        }
        return nearest_node;
    }

    public static Request load_request(int origin, int dest, int submission_t, int[][] dist, int max_wait_t,
                                       float max_tt_p, int max_tt_min, int max_tt_add) {
        int direct_tt = dist[origin][dest];
        int latest_arr_t;
        if (max_tt_p * direct_tt + max_tt_min > max_tt_add) {
            latest_arr_t = (int) (submission_t + max_wait_t + direct_tt + max_tt_add);
        } else {
            latest_arr_t = (int) (submission_t + max_wait_t + (1 + max_tt_p) * direct_tt + max_tt_min);
        }
        Request request = new Request(origin, dest, submission_t, latest_arr_t, null, null,
                "submitted", new ArrayList<>(), false);
        return request;
    }
}
