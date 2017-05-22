package roadgraph;

import java.util.HashSet;
import java.util.Set;

import geography.GeographicPoint;

class MapGraphNode implements Comparable {
	private HashSet<MapGraphEdge> edges;

	private GeographicPoint location;

	private double distance;

	private double actualDistance;

	MapGraphNode(GeographicPoint loc) {
		location = loc;
		edges = new HashSet<MapGraphEdge>();
		distance = 0.0;
		actualDistance = 0.0;
	}

	void addEdge(MapGraphEdge edge) {
		edges.add(edge);
	}

	Set<MapGraphNode> getNeighbors() {
		Set<MapGraphNode> neighbors = new HashSet<MapGraphNode>();
		for (MapGraphEdge edge : edges) {
			neighbors.add(edge.getOtherNode(this));
		}
		return neighbors;
	}

	GeographicPoint getLocation() {
		return location;
	}

	Set<MapGraphEdge> getEdges() {
		return edges;
	}

	public boolean equals(Object o) {
		if (!(o instanceof MapGraphNode) || (o == null)) {
			return false;
		}
		MapGraphNode node = (MapGraphNode) o;
		return node.location.equals(this.location);
	}

	public int HashCode() {
		return location.hashCode();
	}

	public double getDistance() {
		return this.distance;
	}

	public void setDistance(double distance) {
		this.distance = distance;
	}

	public double getActualDistance() {
		return this.actualDistance;
	}

	public void setActualDistance(double actualDistance) {
		this.actualDistance = actualDistance;
	}

	public int compareTo(Object o) {
		MapGraphNode m = (MapGraphNode) o;
		return ((Double) this.getDistance()).compareTo((Double) m.getDistance());
	}

}
