package edu.nyu.intramodal;

import java.util.List;

public class Network {

    int n_nodes;
    int[][] dist;
//    List<Integer>[][] shortest_path;

    public Network(int n_nodes) {
        this.n_nodes = n_nodes;
        //initializing adjacency matrix to the appropriate size
        this.dist = new int[n_nodes][n_nodes];
//        this.shortest_path = new List[n_nodes][n_nodes];
    }

    public void addShortestPath(int src, int dest, int travel_time) {
        dist[src][dest] = travel_time;
//        shortest_path[src][dest] = path;
    }
}
