# ridepooling-with-transfers

This project is a dynamic ride-pooling simulator which allows passengers to transfer between vehicles. It uses a proposed two-deminesional non-myopic cost function approximation policy to make sequential decisions and find the best matches between riders and vehicles in the system. The project is in collaboration with [MOIA](https://www.moia.io), a ride-pooling service operating in Germany. 

![image](https://github.com/BUILTNYU/ridepooling-with-transfers/assets/66441622/34b8a899-b16e-42d7-a809-5ce25cdc99e3)

Ride-pooling services promise to be more efficient than private motorized mobility or ride-hailing. It can reduce costs for both passengers and operators and reduce traffic congestion and environmental impacts. However, achieving high ridership while maintaining system efficiency is still challenging in a ride-pooling service. Allowing transfers within the system is a potential way to improve service availability for customers and increase fleet efficiency. Transfers may enable the system to save on costs and offer service to more users. The static pickup and delivery problem with transfers is an NP-hard problem, while efficient heuristics are necessary for online applications due to the "curse of dimensionality". Myopic decision-making is another well-known issue in dynamic routing, which is often ignored. It is shown that there are opportunity costs to decisions made in routing at a point in time that can impact outcomes in the future horizon. Therefore, it is crucial to adopt a look-ahead approach in fleet operations.

This simulator uses an online policy and algorithm for operating a ride-pooling service with en-route transfers. It can handle large-scale networks with a large number of transfer stops. Two dimensions are considered in the non-myopic policy with transfers: control for opportunity costs due to (1) commitment to serving existing passengers (as observed in conventional dynamic routing) and (2) a new dimension not yet observed in the literature involving a spatiotemporal commitment to meet at a transfer point. The dispatching algorithm decides whether a request should be served by a single vehicle or two. In the latter case, the passenger has one transfer; the passenger is picked up by the first vehicle and dropped off at a transfer stop. The second vehicle picks up the passenger from the transfer stop and drops them off at their destination. The two vehicles, the pickup times, and the place the transfer takes place should be determined by the system. Due to the inconvenience of transfers for passengers, the maximum number of transfers for each passenger is assumed to be one. By adjusting different parameters, the simulator allows implementing and comparing different operations such as a myopic operation without transfers vs a non-myopic operation. 

For more details, please refer to the following paper:
https://www.sciencedirect.com/science/article/pii/S0968090X24001189?via%3Dihub

## Instructions
Simulation parameters and control variables can be determined in the code. Different input and output files are explained below:

### Input files
Example input files are available for the Sioux Falls network. The input files include:
- Network links, which contains the travel time on each link.
- Network nodes, which contains all the nodes in the network and their projected coordinates.
- Service stops, which contains the information of the nodes where passengers can be picked up or dropped off.
- Network travel times, which contains precalculated shortest travel time from each service stop to another.
- Transfer stops, which contains the information of the nodes where transfer activities can take place (optional; the simulator can randomly select a subset of service stops as transfer stops given a ratio between 0 and 1).
- Hub locations, which contains the location of hubs from which vehicles start their operations (optional; the simulator can select the initial location of vehicles randomly from the network nodes).
- Requests, which contains the information of each submitted request including origin, destination, and submission time (optional; the simulator can generate random requests given the total number of submitted requests during operating hours).

### Output files
The output files include:
- Performance evaluation, which contains the values of performance metrics, including average passengers’ in-vehicle time, wait time, perceived journey time (includes in-vehicle time and wait time multiplied by a factor of 1.4), total vehicles’ travel time excluding the dwell times, total vehicles’ empty travel time which refers to the duration that vehicles are moving while not carrying any passengers, average vehicle occupancy, number of rejected requests, number of transfers, simulation runtime.
- Vehicle occupancy information at each time step.
- Requests' information, which includes each request's submission time, pickup time, drop-off time, whether they had transfer or not, vehicles' arrival time at transfer stop in case of having transfers.
- Utilized transfer stops, which contains the information of stops that were used for transfer activities.
