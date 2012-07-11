package statalign.postprocess.utils;

/**
 * A container class that stores a loglikelihood value and a logical value
 * if the sample comes from the burn-in period or not.
 * 
 * @author miklos
 *
 */
public class LogLikelihoodTraceContainer {

	/**
	 * The log-likelihood value
	 */
	public double loglikelihood;
	/**
	 * It tells if the sample comes from the burn-in phase or not.
	 */
	public boolean burnin;
	
	/**
	 * This constructor sets the values in the container.
	 * @param loglikelihood The loglikelihood value
	 * @param burnin It is true, if the sample comes from the burn-in phase.
	 */
	public LogLikelihoodTraceContainer(double loglikelihood, boolean burnin){
		this.loglikelihood = loglikelihood;
		this.burnin = burnin;
		
	}
	
}
