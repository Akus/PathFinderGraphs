package roadgraph;

import geography.GeographicPoint;

class MapGraphEdge {
	/** The name of the road */
	private String roadName;

	/** The type of the road */
	private String roadType;

	/** The two endpoints of the edge */
	private MapGraphNode start;
	private MapGraphNode end;

	/** The length of the road segment, in km */
	private double length;

	static final double DEFAULT_LENGTH = 0.01;

	/**
	 * Create a new MapEdge object
	 * 
	 * @param roadName
	 * @param n1
	 *            The point at one end of the segment
	 * @param n2
	 *            The point at the other end of the segment
	 * 
	 */
	MapGraphEdge(String roadName, MapGraphNode n1, MapGraphNode n2) {
		this(roadName, "", n1, n2, DEFAULT_LENGTH);
	}

	MapGraphEdge(String roadName, String roadType, MapGraphNode n1, MapGraphNode n2) {
		this(roadName, roadType, n1, n2, DEFAULT_LENGTH);
	}

	MapGraphEdge(String roadName, String roadType, MapGraphNode n1, MapGraphNode n2, double length) {
		this.roadName = roadName;
		start = n1;
		end = n2;
		this.roadType = roadType;
		this.length = length;
	}

	String getEdgeName() {
		return roadName;
	}
	
	MapGraphNode getEndNode() {
		return end;
	}

	GeographicPoint getStartPoint() {
		return start.getLocation();
	}

	GeographicPoint getEndPoint() {
		return end.getLocation();
	}

	double getLength() {
		return length;
	}

	public String getRoadName() {
		return roadName;
	}

	MapGraphNode getOtherNode(MapGraphNode node) {
		if (node.equals(start))
			return end;
		else if (node.equals(end))
			return start;
		throw new IllegalArgumentException("Looking for " + "a point that is not in the edge");
	}

}