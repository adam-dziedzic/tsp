import static org.apache.commons.math3.util.ArithmeticUtils.binomialCoefficient;

import java.io.File;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Tsp2 {

	/** Logging setup. */
	private Logger LOGGER = Logger.getLogger(tsp.class.getName());

	/** Level of the logging. */
	private Level LOGGING_LEVEL = Level.OFF;

	/** Tsp needs only previous subproblem to compute current subproblem. */
	private double[][] previousSubproblem;

	/** Tsp computes currentSubproblem from the previous subproblem. */
	private double[][] currentSubproblem;

	/** Name of the file with data. */
	private String fileName;

	/** Collection of the cities (nodes). */
	private Node[] nodes;

	/** Number of nodes (cities) for the tsp problem. */
	private int nodesSize;

	/** Distances between nodes. */
	private double[][] distances;

	/**
	 * Set nodes size.
	 * 
	 * @param nodesSize
	 *            - number of nodes
	 */
	public void setNodesSize(int nodesSize) {
		this.nodesSize = nodesSize;
	}

	/**
	 * start the program
	 * 
	 * @param args
	 *            [0] name of the file with data
	 */
	public static void main(String[] args) {
		String fileName = "tsp4.txt";
		if (args.length > 0) {
			fileName = args[0];
		}
		Tsp2 tsp = new Tsp2(fileName);
		tsp.start();
	}

	/**
	 * Create tsp object to solve the problem.
	 * 
	 * @param fileName
	 *            - name of the file with data.
	 */
	public Tsp2(String fileName) {
		this.fileName = fileName;
	}

	/**
	 * Start the program to compute the shortest graph cycle in tsp problem.
	 * 
	 */
	public void start() {
		LOGGER.setLevel(LOGGING_LEVEL);
		readNodes();
		initDistances();

		long startTime = System.currentTimeMillis();
		initCurrentSubproblem();
		double result = tspAlgorithm();
		long endTime = System.currentTimeMillis();

		System.out.println("Result: " + result);
		System.out.println("Elapsed time: " + (endTime - startTime) / 1000.0
		        + " seconds.");
	}

	/**
	 * TSP dynamic algorithm.
	 * 
	 * Computes only the shortest path without giving the path itself. This is
	 * why the final complexity is O(n^2*2^n) instead of O(n!);
	 */
	public double tspAlgorithm() {
		// setSize corresponds to m in the original algorithm
		for (int setSize = 2; setSize < nodesSize; ++setSize) {
			previousSubproblem = currentSubproblem;
			// LOGGER.info(show2Darray(previousSubproblem));
			currentSubproblem = new double[(int) binomialCoefficient(
			        nodesSize - 1, setSize)][];

			// init value for Gosper's generator
			int Sint = getGosperNumber(nodesSize - 1, setSize, 0);
			while (Sint != -1) {
				double[] distancesS = new double[setSize];
				int[] S = getSet(Sint, setSize);
				for (int i = 0; i < S.length; ++i) {
					int j = S[i];
					LOGGER.info("J: " + j);
					int[] Swithoutj = arrayWithoutIndexI(S, i);
					double[] distancesSWithoutj = previousSubproblem[getIndex(Swithoutj)];
					double minValue = Double.MAX_VALUE;
					for (int index = 0; index < distancesSWithoutj.length; ++index) {
						int k = Swithoutj[index];
						double distanceFromEndToK = distancesSWithoutj[index];
						double distanceFromKToJ = distances[k][j];
						double distanceFromEndToJ = distanceFromEndToK
						        + distanceFromKToJ;
						if (distanceFromEndToJ < minValue) {
							minValue = distanceFromEndToJ;
						}
					}
					distancesS[i] = minValue;
				}
				currentSubproblem[getIndex(S)] = distancesS;
				Sint = getGosperNumber(nodesSize - 1, setSize, Sint);
			}
			// LOGGER.info(show2Darray(currentSubproblem));
		}
		return getMinValue();
	}

	/**
	 * Get final min Value of the tsp dynamic programming algorithm after
	 * closing the cycle;
	 * 
	 * @return minimal value of the path for the salesman
	 */
	public double getMinValue() {
		double minGlobal = Double.MAX_VALUE;
		int Sint = getGosperNumber(nodesSize - 1, nodesSize - 1, 0);
		int[] S = getSet(Sint, nodesSize - 1);
		for (int i = 0; i < currentSubproblem[0].length; ++i) {
			double finalDistance = currentSubproblem[0][i]
			        + distances[S[i]][nodesSize - 1];
			// System.out.println("Final distance: " + S[i]);
			if (finalDistance < minGlobal)
				minGlobal = finalDistance;
		}
		// LOGGER.info(show2Darray(currentSubproblem));
		return minGlobal;
	}

	/**
	 * Remove element with index i from array S.
	 * 
	 * @param S
	 *            array of int elements
	 * @param j
	 *            index of the element to be removed from array S
	 */
	public int[] arrayWithoutIndexI(int[] S, int i) {
		int[] result = new int[S.length - 1];
		int index = 0; // index of the new array result to insert elements
		for (int j = 0; j < i; ++j) {
			result[index++] = S[j];
		}
		for (int j = i + 1; j < S.length; ++j) {
			result[index++] = S[j];
		}
		return result;
	}

	public void initCurrentSubproblem() {
		currentSubproblem = new double[(int) binomialCoefficient(nodesSize - 1,
		        1)][];
		// initial subset for the subproblem of size 1
		int Sint = getGosperNumber(nodesSize - 1, 1, 0);
		while (Sint != -1) {
			int[] S = getSet(Sint, 1);
			// distances list for the first set
			double[] distancesList = new double[1];
			// distance to the last node;
			// S[0] - first element in the initial set
			distancesList[0] = distances[S[0]][nodesSize - 1];
			LOGGER.info("Get index S: " + getIndex(S));
			currentSubproblem[getIndex(S)] = distancesList;
			// generate the next Gosper's number
			Sint = getGosperNumber(nodesSize - 1, 1, Sint);
		}
	}

	/**
	 * Gosper's hack to generate subsets. Generate subsets of a set of size n,
	 * each subset of size k
	 * 
	 * We use only (size - 1) bits because "0" (zero) represents node 0 which is
	 * always included in the set S.
	 * 
	 * @return Goseper's number from k combinations
	 */
	public int getGosperNumber(int n, int k, int last) {
		int set = last;
		if (last == 0) { // initialize Gosper's generator
			set = (1 << k) - 1;
			return set;
		}
		int limit = (1 << n);
		// Gosper's hack:
		int c = set & -set;
		int r = set + c;
		set = (((r ^ set) >>> 2) / c) | r;

		if (set < limit) {
			return set;
		}
		return -1;
	}

	/**
	 * Initialize distances between nodes and put them in the distances table.
	 */
	public void initDistances() {
		distances = new double[nodesSize][nodesSize];
		for (int i = 0; i < nodesSize; ++i) {
			for (int j = i + 1; j < nodesSize; ++j) {
				double distance = nodes[i].getDistance(nodes[j]);
				distances[i][j] = distance;
				distances[j][i] = distance;
			}
		}
	}

	/**
	 * Format of the file: The first line indicates the number of cities. Each
	 * city is a point in the plane, and each subsequent line indicates the x-
	 * and y-coordinates of a single city.
	 * 
	 * Example of a file: 25 20833.3333 17100.0000 20900.0000 17066.6667
	 */
	public void readNodes() {
		Scanner sc = null;
		try {
			sc = new Scanner(new File(fileName));
		} catch (Exception e) {
			e.printStackTrace();
		}
		nodesSize = sc.nextInt();
		nodes = new Node[nodesSize];
		int i = 0;
		while (sc.hasNextDouble()) {
			double x = sc.nextDouble();
			double y = sc.nextDouble();
			Node node = new Node(x, y);
			nodes[i++] = node;
		}
	}

	/**
	 * Get a set of numbers from the given int n.
	 * 
	 * @param gosperNumber
	 *            - int number with set bits
	 * @param size
	 *            - number of significant bits
	 * 
	 * @return - set of number corresponding to "1" bits in the integer value
	 */
	public int[] getSet(int gosperNumber, int setSize) {
		int[] set = new int[setSize];
		int index = 0;
		// number of nodes without last node
		// nodes.size() - 1 == number of significant bits
		int size = nodesSize - 1;
		for (int j = size - 1; j >= 0; --j) {
			if (((gosperNumber & 1 << j) != 0) {
				set[index++] = j; // bits: 101 -> nodes: {2,0}
			}
		}
		return set;
	}

	/**
	 * For a given set of ints return Gosper's int.
	 * 
	 * @param set
	 *            set of int numbers
	 * @return gosper's int
	 */
	public int getInt(int[] set) {
		int result = 0;
		for (int i : set) {
			result += (1 << i);
		}
		return result;
	}

	/**
	 * Function mapping combinations to natural numbers.
	 * 
	 * Numbering of the k-combinations of a given set. That is, you want a
	 * function whose output is equivalent to that of the following. Given a set
	 * A of size n, and some number k:
	 * 
	 * enumerate all subsets of A of size k (there are binomial(n,k) of them),
	 * Arrange them in some (any) order, Give them numbers from 0 to
	 * binomial(n,k)−1.
	 * 
	 * It is possible to efficiently compute such a function, and the most
	 * common one is known as "combinatorial number system". The function is as
	 * follows:
	 * 
	 * Label the elements of A as integers 0 to n−1. Given a subset S of A of
	 * size k, write the elements of the subset in decreasing order, as
	 * S={ck,ck−1,…,c1} with ck>⋯>c1≥0. The function
	 * f(S)=binomial(ck,k)+⋯+binomial(c2,2)+binomial(c1,1).
	 * 
	 * This is quite fast to compute, and it's also fast to go in the other
	 * direction (find the set S given f(S)).
	 * 
	 * There is a similar factorial numbering system for numbering permutations.
	 * 
	 * source: http://bit.ly/15eLoQQ
	 * 
	 * A = {0,1, ... , nodesSize}
	 * 
	 * @param set
	 *            - sorted array S (decreasing order)
	 * @return index of the set (combination)
	 */
	public int getIndex(int[] set) {
		int result = 0;
		int k = set.length; // subset S of A of size k
		for (int c : set) {
			if (c >= k) {
				result += binomialCoefficient(c, k--);
			}
		}
		return result;
	}

	/**
	 * Show two dimensional array.
	 * 
	 * @param tab
	 *            two dimensional array
	 */
	public String show2Darray(double[][] tab) {
		StringBuilder s = new StringBuilder();
		for (double[] row : tab) {
			for (double n : row) {
				s.append(n + " ");
			}
			s.append("\n");
		}
		return s.toString();
	}

}
