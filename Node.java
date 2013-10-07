/**
 * Node representing a city on a map.
 * 
 */
class Node {
	private double x;
	private double y;

	public Node(double x, double y) {
		this.x = x;
		this.y = y;
	}

	public double getDistance(Node n) {
		double xSquared = Math.pow(x - n.x, 2);
		double ySquared = Math.pow(y - n.y, 2);
		return Math.sqrt(xSquared + ySquared);
	}

	public String toString() {
		return "(x: " + x + " y: " + y + ")";
	}
}