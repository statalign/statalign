package statalign.postprocess.utils;

/**
 * A container class that stores a parameter values from the StructAlign plugin
 * and a logical value indicating whether the sample comes from the burn-in period or not.
 * 
 * @author herman
 *
 */
public class StructAlignTraceParameters {

	/**
	 * The log-likelihood value
	 */
	public double[] sigma2;
	public double sigma2Hier;
	public double nu;
	public double tau;
	public boolean globalSigma = true;
	public double epsilon;
	/**
	 * It tells if the sample comes from the burn-in phase or not.
	 */
	public boolean burnin;
	
	/**
	 * This constructor sets the values in the container.
	 * @param loglikelihood The loglikelihood value
	 * @param burnin It is true, if the sample comes from the burn-in phase.
	 */
	public StructAlignTraceParameters(boolean burnin){
//		this.loglikelihood = loglikelihood;
		this.burnin = burnin;
//		
	}
}
