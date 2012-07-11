package statalign.base.hmm;

/**
 * 
 * This is the abstract class of pair-HMMs that emit characters into one unobservable
 * sequence, namely that aligns two observable sequences via an unobservable
 * ancestral sequence.
 * 
 * @author novak
 *
 */
abstract public class HmmSilent extends Hmm {
	
	/**
	 * Calculates the transition matrix given edge lengths t1 and t2 of a cherry-motif.
	 * 
	 * @param transMatrix If null, the function allocates memory for the return matrix
	 *        If not, this is used for storing the matrix. For speeding-up purposes,
	 *        note that allocating memory is time-consuming. In this way we need to allocate
	 *        memory when asked.
	 * @param t1 The length of the left edge
	 * @param t2 The length of the right edge.
	 * @return The calculated transition matrix.
	 */
	public abstract double[][] preCalcTransMatrix(double[][] transMatrix, double t1, double t2);

	/**
	 * Silent state elimination.
	 * 
	 * Transforms transMatrix into the transition matrix obtained by eliminating the
	 * silent state of transMatrix.
	 * 
	 * @param redTransMatrix Array to store the results, may be left null
	 * @param transMatrix Input transition matrix with a silent state
	 * @return The calculated ('reduced') transition matrix (redTransMatrix or a newly
	 *   allocated array if it is null)
	 */
	public abstract double[][] preCalcRedTransMatrix(double redTransMatrix[][], double transMatrix[][]);
	
	/**
	 * Returns the index of the silent state.
	 * 
	 * @return 	Index of the silent state
	 */
	public abstract int getSilent();
	
}
