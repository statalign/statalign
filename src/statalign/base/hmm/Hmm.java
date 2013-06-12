package statalign.base.hmm;

/**
 * The abstract class of Hidden Markov Models.
 * 
 * @author novak
 *
 */
public abstract class Hmm {

	/**
	 * HMM model parameters. Implementing classes can use their own assignment and
	 * any number of parameters.
	 * In the case of the TKF92 model (HmmTkf92), values are: r, lambda and mu.
	 */
	public double[] params;
	
   /**
    * Abstract function that specifies the emission of each state in the HMM in the form
    * of an array A: A[s][i]=1 iff state i emits into sequence s (s=0 is the parent sequence;
    * for pair-HMMs s=1 is the child sequence, for 3-seq HMMs s=1 is left and s=2 is right
    * child).
    * HMM implementations must override this and return their own emission patterns.
    * @return array specifying emission patterns as described above
    */
	public abstract int[][] getStateEmit();
	
	/**
	 * Abstract function returning an array A that helps identify the state with some
	 * emission pattern. The emission pattern is coded as a single integer: e.g. for a 3-seq
	 * HMM, emission into the parent sequence and right child only is coded as 101 binary,
	 * i.e. 4+1=5. Then, A[5] gives the index of the state with the above emission pattern.
	 * HMM implementations must override this and return an array corresponding to their own
	 * state - emission pattern assignment.
	 * @return array describing the emission pattern to state conversion
	 */
	public abstract int[] getEmitPatt2State();
	
	/**
	 * Returns the index of the start state.
	 * 
	 * @return The index of the start state;
	 */
	public abstract int getStart();
	
	/**
	 * Returns the index of the end state.
	 * 
	 * @return The index of the end state;
	 */
	public abstract int getEnd();
	
	/**
	 * Returns the logarithm of the stationary probability of generating
	 * a sequence of <code>length</code> characters under the HMM.
	 * @param length The length of the sequence whose stationary
	 * probability is to be computed.
	 * @return The logarithm of the stationary probability under the HMM.
	 */
	public double getLogStationaryProb(int length) { return 0.0; }
	
}
