/**
 * @author UCSD MOOC development team and YOU
 * 
 * A class which reprsents a graph of geographic locations
 * Nodes in the graph are intersections between 
 *
 */
package roadgraph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Consumer;

import geography.GeographicPoint;
import util.GraphLoader;

/**
 * @author UCSD MOOC development team
 *
 *         A class which represents a graph of geographic locations Nodes in the
 *         graph are intersections of multiple roads. Edges are the roads.
 *
 */
public class MapGraph {

	private HashMap<GeographicPoint, MapGraphNode> geoNodeMap;
	private HashSet<MapGraphEdge> graphEdges;

	/**
	 * Create a new empty MapGraph
	 *
	 */
	public MapGraph() {
		geoNodeMap = new HashMap<GeographicPoint, MapGraphNode>();
		graphEdges = new HashSet<MapGraphEdge>();
	}

	/**
	 * Get the number of vertices (road intersections) in the graph
	 * 
	 * @return The number of vertices in the graph.
	 */
	public int getNumVertices() {
		return geoNodeMap.values().size();
	}

	/**
	 * Get the number of road segments in the graph
	 * 
	 * @return The number of edges in the graph.
	 */
	public int getNumEdges() {
		return graphEdges.size();
	}

	/**
	 * Add a node corresponding to an intersection
	 *
	 * @param latitude
	 *            The latitude of the location
	 * @param longitude
	 *            The longitude of the location
	 */
	public void addVertex(double latitude, double longitude) {
		GeographicPoint pt = new GeographicPoint(latitude, longitude);
		this.addVertex(pt);
		System.out.println("The vertex: " + pt + " has been successfully added to the graph. It's current value is: "
				+ geoNodeMap.get(pt));
	}

	/**
	 * Add a node corresponding to an intersection at a Geographic Point
	 *
	 * @param location
	 *            The location of the intersection
	 */
	public void addVertex(GeographicPoint location) {
		MapGraphNode node = geoNodeMap.get(location);
		if (node == null) {
			node = new MapGraphNode(location);
			geoNodeMap.put(location, node);
		} else {
			System.out.println("This location is already occupied by a Node: " + location);
		}

	}

	// removing a Vertex made by Akos
	public void removeVertex(double latitude, double longitude) {

		GeographicPoint location = new GeographicPoint(latitude, longitude);

		MapGraphNode node = geoNodeMap.get(location);
		if (node != null) {
			node = null;
			geoNodeMap.put(location, node);
			System.out.println("This location has been removed successfully: " + location
					+ " The current value of the node is now: " + node);

		} else {
			System.out.println("This location is already empty: " + location);
		}

	}

	/**
	 * Add an edge representing a segment of a road. Precondition: The
	 * corresponding Nodes must have already been added to the graph.
	 * 
	 * @param roadName
	 *            The name of the road
	 * @param roadType
	 *            The type of the road
	 */
	public void addEdge(double lat1, double lon1, double lat2, double lon2, String roadName, String roadType) {
		// Find the two Nodes associated with this edge.
		GeographicPoint from = new GeographicPoint(lat1, lon1);
		GeographicPoint to = new GeographicPoint(lat2, lon2);

		MapGraphNode fromNode = geoNodeMap.get(from);
		MapGraphNode toNode = geoNodeMap.get(to);

		// validation
		if (fromNode == null)
			throw new NullPointerException(from + " is not in graph");
		if (toNode == null)
			throw new NullPointerException(to + " is not in graph");

		addEdge(fromNode, toNode, roadName, roadType, MapGraphEdge.DEFAULT_LENGTH);

	}

	public void addEdge(GeographicPoint from, GeographicPoint to, String roadName, String roadType) {

		MapGraphNode fromNode = geoNodeMap.get(from);
		MapGraphNode toNode = geoNodeMap.get(to);

		if (fromNode == null)
			throw new NullPointerException(from + "is not in graph");
		if (toNode == null)
			throw new NullPointerException(to + "is not in graph");

		addEdge(fromNode, toNode, roadName, roadType, MapGraphEdge.DEFAULT_LENGTH);
	}

	public void addEdge(GeographicPoint from, GeographicPoint to, String roadName, String roadType, double length) {
		MapGraphNode fromNode = geoNodeMap.get(from);
		MapGraphNode toNode = geoNodeMap.get(to);

		if (fromNode == null)
			throw new NullPointerException(from + "is not in graph");
		if (toNode == null)
			throw new NullPointerException(to + "is not in graph");

		addEdge(fromNode, toNode, roadName, roadType, length);
	}

	public boolean nodeChecker(GeographicPoint point) {
		return geoNodeMap.containsKey(point);
	}

	private void addEdge(MapGraphNode fromNode, MapGraphNode toNode, String roadName, String roadType, double length) {
		MapGraphEdge myEdge = new MapGraphEdge(roadName, roadType, fromNode, toNode, length);
		graphEdges.add(myEdge);
		fromNode.addEdge(myEdge);
	}

	public Collection<GeographicPoint> getVertices() {
		return geoNodeMap.keySet();
	}

	private Set<MapGraphNode> getNeighbors(MapGraphNode node) {
		return node.getNeighbors();
	}

	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		Consumer<GeographicPoint> temp = (x) -> {
		};
		return bfs(start, goal, temp);
	}

	/**
	 * Find the path from start to goal using Breadth First Search
	 *
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched) {

		if (!validatorHelper(start, goal)) {
			return null;
		}

		MapGraphNode startNode = geoNodeMap.get(start);
		MapGraphNode endNode = geoNodeMap.get(goal);

		HashMap<MapGraphNode, MapGraphNode> parentMap = new HashMap<MapGraphNode, MapGraphNode>();
		Queue<MapGraphNode> toExplore = new LinkedList<MapGraphNode>();
		HashSet<MapGraphNode> visited = new HashSet<MapGraphNode>();
		toExplore.add(startNode);
		MapGraphNode next = null;

		while (!toExplore.isEmpty()) {
			next = toExplore.remove();

			// hook for visualization
			nodeSearched.accept(next.getLocation());

			if (next.equals(endNode))
				break;

			for (MapGraphNode neighbor : getNeighbors(next)) {
				if (!visited.contains(neighbor)) {
					visited.add(neighbor);
					parentMap.put(neighbor, next);
					toExplore.add(neighbor);
				}
			}
		}

		return reconstructPath(parentMap, startNode, endNode, next.equals(endNode));
	}

	private boolean validatorHelper(GeographicPoint start, GeographicPoint goal) {
		if (start == null || goal == null) {
			throw new NullPointerException("Null node error!");
		}
		if (geoNodeMap.get(start) == null) {
			System.err.println(start + " not exist");
			return false;
		}
		if (geoNodeMap.get(goal) == null) {
			System.err.println(goal + " not exist");
			return false;
		}
		return true;
	}

	private List<GeographicPoint> reconstructPath(HashMap<MapGraphNode, MapGraphNode> parentMap, MapGraphNode start,
			MapGraphNode goal, boolean pathFound) {
		if (!pathFound) {
			System.out.println("No path from " + start + " to " + goal);
			return null;
		}
		LinkedList<GeographicPoint> path = new LinkedList<GeographicPoint>();
		MapGraphNode current = goal;

		while (!current.equals(start)) {
			path.addFirst(current.getLocation());
			current = parentMap.get(current);
		}

		path.addFirst(start.getLocation());
		return path;
	}

	/**
	 * Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		// You do not need to change this method.
		Consumer<GeographicPoint> temp = (x) -> {
		};
		return dijkstra(start, goal, temp);
	}

	/**
	 * Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @param nodeSearched
	 *            A hook for visualization. See assignment instructions for how
	 *            to use it.
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched) {
		return ultimateSearch(start, goal, nodeSearched, (a, b) -> 0.0);
	}

	private List<GeographicPoint> ultimateSearch(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched, BiFunction<MapGraphNode, MapGraphNode, Double> f) {

		if (!validatorHelper(start, goal)) {
			return null;
		}

		MapGraphNode startNode = geoNodeMap.get(start);
		MapGraphNode endNode = geoNodeMap.get(goal);

		HashMap<MapGraphNode, MapGraphNode> parentMap = new HashMap<MapGraphNode, MapGraphNode>();
		PriorityQueue<MapGraphNode> toExplore = new PriorityQueue<MapGraphNode>();
		HashSet<MapGraphNode> visited = new HashSet<MapGraphNode>();

		initializeDistances();

		startNode.setActualDistance(0.0);
		startNode.setDistance(0.0);

		toExplore.add(startNode);
		MapGraphNode next = null;

		while (!toExplore.isEmpty()) {
			next = toExplore.poll();

			if (!visited.contains(next)) {
				visited.add(next);

				// hook for visualization
				nodeSearched.accept(next.getLocation());

				if (next.equals(endNode)) {
					break;
				}

				HashMap<MapGraphNode, Double> distancesMap = calculateDistanesMap(next);

				for (MapGraphNode neighbor : getNeighbors(next)) {
					if (!visited.contains(neighbor)) {
						double distanceOfNode = next.getActualDistance() + distancesMap.get(neighbor);
						if (distanceOfNode < neighbor.getActualDistance()) {
							neighbor.setActualDistance(distanceOfNode);
							distanceOfNode += f.apply(neighbor, endNode);
							neighbor.setDistance(distanceOfNode);
							parentMap.put(neighbor, next);
							toExplore.offer(neighbor);
						}
					}
				}
			}
		}

		System.out.println("Visited: " + visited.size());

		// Reconstruct the parent path
		return reconstructPath(parentMap, startNode, endNode, endNode.equals(next));
	}

	private void initializeDistances() {
		for (MapGraphNode m : geoNodeMap.values()) {
			m.setActualDistance(Double.MAX_VALUE);
			m.setDistance(Double.MAX_VALUE);
		}
	}

	private HashMap<MapGraphNode, Double> calculateDistanesMap(MapGraphNode next) {
		HashMap<MapGraphNode, Double> distancesMap = new HashMap<>();

		for (MapGraphEdge e : next.getEdges()) {
			distancesMap.put(e.getEndNode(), e.getLength());
		}
		return distancesMap;
	}

	/**
	 * Find the path from start to goal using A-Star search
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		Consumer<GeographicPoint> temp = (x) -> {
		};
		return aStarSearch(start, goal, temp);
	}

	/**
	 * Find the path from start to goal using A-Star search
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @param nodeSearched
	 *            A hook for visualization. See assignment instructions for how
	 *            to use it.
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched) {
		return ultimateSearch(start, goal, nodeSearched, (a, b) -> a.getLocation().distance(b.getLocation()));
	}

	// my Awesome Greedy TSP

	public List<GeographicPoint> greedyTSPtest(GeographicPoint homeTownPoint) {

		System.out.println("Location of homeTown node: " + geoNodeMap.get(homeTownPoint).getLocation());

		MapGraphNode homeTownNode = geoNodeMap.get(homeTownPoint);

		// current = Hometown
		MapGraphNode currentNode = homeTownNode;
		GeographicPoint currentNodePoint = homeTownPoint;
		// nodesToVisit = all other nodes

		// adding all nodes to nodesToVisit

		Collection<GeographicPoint> allNodePoints = this.getVertices();

		System.out.println(allNodePoints);

		// removing safely HomeTownPoint from Collection with Iterator
		Iterator<GeographicPoint> setIterator = allNodePoints.iterator();

		while (setIterator.hasNext()) {

			GeographicPoint currentElement = setIterator.next();
			if (currentElement.equals(homeTownPoint)) {
				setIterator.remove();
			}
		}

		Collection<GeographicPoint> nodesToVisit = allNodePoints;
		System.out.println("nodesToVisit after removing homeTownPoint: " + nodesToVisit);

		// while (more nodes to visit) -> O(n-1)
		// selected_node = closest to current -> O(n)
		// add selected_node to bestPath
		// remove current node from nodes_to_visit
		// current = selected node -> total => O(n-1)*n -> Order of n squared

		List<GeographicPoint> bestPath = new ArrayList();
		MapGraphNode closestNode = null;
		GeographicPoint closestNodePoint = null;
		System.out.println("nodesToVisit size: " + nodesToVisit.size() + "\nnodesToVisit isEmpty? " + nodesToVisit.isEmpty());
		
		while (nodesToVisit.size() != 0) {

			// finding the shortest edge

			double shortestDist = Double.POSITIVE_INFINITY;
			
			
			System.out.println("CurrentNodePoint before FOR CYCLE: " + currentNodePoint);
			System.out.println("Current Node: location before FOR CYCLE: " + currentNode.getLocation());
			System.out.println("Edges of currentNode: " + currentNode.getEdges());
			
			for (MapGraphEdge edge : currentNode.getEdges()) {

				System.out.println("\n###########################");
				System.out.println("CurrentNodePoint: " + currentNodePoint);
				System.out.println("Current Node: location: " + currentNode.getLocation());
				
				
				System.out.println("Current Edge: name: " + edge.getEdgeName());
				System.out.println("Current Edge: length: " + edge.getLength());

				if (edge.getLength() < shortestDist) {

					shortestDist = edge.getLength();
					System.out.println("shortestDist now is: " + shortestDist);
					closestNode = edge.getEndNode();
					closestNodePoint = closestNode.getLocation();
					System.out.println("This edge is now the shortest: " + edge.getEdgeName() + " It's length is: "
							+ edge.getLength());
					System.out.println("The closest node to " + currentNode.toString() + " is now: " + closestNode);
					System.out.println("The GeographicPoints of closest node is now: " + closestNodePoint);

				}
				System.out.println("End of inner for cycle");
				System.out.println("###########################\n");

			}

			System.out.println("Line# 490: closestNode is: " + closestNode);
			System.out.println("closestNodePoint is: " + closestNodePoint);

			// adding closestNodePoint to bestPath
			bestPath.add(closestNodePoint);

			System.out.println("bestPath: " + bestPath.toString());
			
			System.out.println("currentNode before closestNode is: " + currentNode);

			currentNode = closestNode;
			currentNodePoint = currentNode.getLocation();

			System.out.println("currentNode after closestNode is: " + currentNode);
			System.out.println("closestNode is: " + closestNode);
			System.out.println("closestNodePoint is: " + closestNodePoint);
			
			System.out.println("currentNodePoint before closestNodePoint is: " + currentNodePoint);
			//currentNodePoint = closestNodePoint;
			System.out.println("currentNodePoint after closestNodePoint is: " + closestNodePoint);

			// removing closestNodePoint from nodesToVisit
			System.out.println("nodesToVisit size: " + nodesToVisit.size() + "\nnodesToVisit isEmpty? " + nodesToVisit.isEmpty());

			Iterator<GeographicPoint> setIterator3 = nodesToVisit.iterator();

			while (setIterator3.hasNext()) {

				GeographicPoint currentElement = setIterator3.next();
				if (currentElement.equals(closestNodePoint)) {
					setIterator3.remove();
				}
			}
			System.out.println("nodesToVisit size: " + nodesToVisit.size() + "\nnodesToVisit isEmpty? " + nodesToVisit.isEmpty());

			System.out.println("nodesToVisit after removing closestNodePoint: " + nodesToVisit);
			System.out.println("End of outer while cycle");
			System.out.println("###########################\n");

		}

		return bestPath;

	}

	public List<MapGraphNode> greedyTSP(double latitude, double longitude) {

		// bestPath = []
		List<MapGraphNode> bestPath = new ArrayList();

		GeographicPoint homeTownPoint = new GeographicPoint(latitude, longitude);

		MapGraphNode homeTownNode = geoNodeMap.get(homeTownPoint);

		// current = Hometown
		MapGraphNode currentNode = homeTownNode;

		// nodesToVisit = all other nodes

		// adding all nodes to nodesToVisit
		Collection<GeographicPoint> nodesToVisit = this.getVertices();

		// removing Hometown from nodesToVisit
		nodesToVisit.remove(homeTownPoint);

		// testing nodesToVisit
		// return nodesToVisit.toString();

		// while (more nodes to visit) -> O(n-1)
		// selected_node = closest to current -> O(n)
		// add selected_node to bestPath
		// remove current node from nodes_to_visit
		// current = selected node -> total => O(n-1)*n -> Order of n squared

		while (!nodesToVisit.isEmpty()) {

			// finding the shortest edge
			// List NeighborEdgeDistances = new ArrayList();

			double shortestDist = Double.POSITIVE_INFINITY;
			MapGraphNode closestNode = null;

			for (MapGraphEdge edge : currentNode.getEdges()) {

				System.out.println("Current Edge: " + edge.getEdgeName());

				if (edge.getLength() < shortestDist) {
					shortestDist = edge.getLength();
					closestNode = edge.getEndNode();
					System.out.println("This edge is now the shortest: " + edge.getEdgeName() + " It's length is: "
							+ (double) edge.getLength());
					System.out.println("The closest node to " + currentNode.toString() + " is now: " + closestNode);

				}

			}
			bestPath.add(closestNode);

			nodesToVisit.remove(closestNode);

			currentNode = closestNode;

		}

		return bestPath;

	}

	// main method for testing
	public static void main(String[] args) {

		// creating the graph manually

		MapGraph simpleGraph = new MapGraph();

		// creating GeographicPoints of nodes
		GeographicPoint pointA = new GeographicPoint(7.0, 1.0);
		GeographicPoint pointB = new GeographicPoint(6.0, 6.0);
		GeographicPoint pointC = new GeographicPoint(1.0, 6.0);
		GeographicPoint pointD = new GeographicPoint(1.0, 1.0);

		// adding nodes
		simpleGraph.addVertex(pointA); // adding A
		simpleGraph.addVertex(pointB); // adding B
		simpleGraph.addVertex(pointC); // adding C
		simpleGraph.addVertex(pointD); // adding D
		
		// WARNING: all nodes must have edges!!! unless it will skip the inner for loop in greedyTSP function!!!
		simpleGraph.addEdge(pointA, pointB, "AB", "simple", 5.0);
		simpleGraph.addEdge(pointB, pointA, "BA", "simple", 5.0);
		simpleGraph.addEdge(pointA, pointC, "AC", "simple", 6.0);
		simpleGraph.addEdge(pointC, pointA, "CA", "simple", 6.0);
		simpleGraph.addEdge(pointA, pointD, "AD", "simple", 6.0);
		simpleGraph.addEdge(pointD, pointA, "DA", "simple", 6.0);
		simpleGraph.addEdge(pointB, pointC, "BC", "simple", 5.0);
		simpleGraph.addEdge(pointC, pointB, "CB", "simple", 5.0);
		simpleGraph.addEdge(pointB, pointD, "BD", "simple", 6.0);
		simpleGraph.addEdge(pointD, pointB, "DB", "simple", 6.0);
		simpleGraph.addEdge(pointD, pointC, "DC", "simple", 5.0);
		simpleGraph.addEdge(pointC, pointD, "CD", "simple", 5.0);

		/*
		 * I have implemented the Greedy Traveling Salesman algorithm. This algorithm starts from the start node and visits all the other nodes only once. Finally it gets back to the start node. Greedy TSP does this job with Order of n squared performance which is not the best but can be optimal.
		 * 
		 * First I tried to implement the algorithm based on the map files and the Loader class but it uses a different distance function which made my job harder.
So I created a simple Graph with my own MapGraphNode class instead.
The second problem was that I missed to link the nodes in both directions so the algorithm couldn't find its way. But I realized that mistake and solved it.
		 */
		
		
		System.out.println("simpleGraph with disjkstra: " + simpleGraph.dijkstra(pointA, pointC));
		System.out.println("simpleGraph with bfs: " + simpleGraph.bfs(pointA, pointC));
		System.out.println("simpleGraph with aStar: " + simpleGraph.aStarSearch(pointA, pointC));
		// simpleGraph.greedyTSP(7.0, 1.0);

		System.out.println("Number of Edges: " + simpleGraph.getNumEdges());
		System.out.println("Nodes: " + simpleGraph.getVertices());

		for (MapGraphEdge e : simpleGraph.graphEdges) {
			System.out.println(e.getEdgeName());
		}

		System.out.println(simpleGraph.getVertices());

		/*
		 * // testing remove Vertex
		 * 
		 * double latitude = 6.0; double longitude = 6.0;
		 * 
		 * Greedy_TSP_map.removeVertex(latitude, longitude);
		 * 
		 * // readding Vertex Greedy_TSP_map.addVertex(latitude, longitude);
		 */

		// calling greedyTSP
		// setting coordinates of Hometown

		GeographicPoint homeTownPoint = new GeographicPoint(7.0, 1.0);

		// testing nodesToVisit
		// System.out.println("The nodesToVisit without Hometown: " +
		// simpleGraph.greedyTSP(latitude, longitude));
		simpleGraph.greedyTSPtest(homeTownPoint);

	}

}