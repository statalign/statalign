package statalign.base.hmm;

/**
 * This is an abstract class for HMMs that are pair-HMMs and do not emit characters into
 * unobservable sequences.
 * 
 * @author novak
 *
 */
public abstract class Hmm2 extends Hmm {
	
	/**
	 * Calculates the transition matrix given an edge length t.
	 * 
	 * @param transMatrix If null, the function allocates memory for the return matrix
	 *        If not, this is used for storing the matrix. For speeding-up purposes,
	 *        note that allocating memory is time-consuming. In this way we need to allocate
	 *        memory when asked.
	 * @param t The edge length parameter.
	 * @return The calculated transition matrix.
	 */
	public abstract double[][] preCalcTransMatrix(double[][] transMatrix, double t);
	
}
