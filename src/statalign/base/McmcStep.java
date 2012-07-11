package statalign.base;

/**
 * Provides information on an MCMC step and additional MCMC statistics.
 * 
 * @author novak
 */
public class McmcStep {

	/** <code>true</code> if this step was proposed in the burn-in period */
	public boolean burnIn;
	
	/** loglikelihood before this step */
	public double oldLogLike;
	/** loglikelihood of proposed step */
	public double proposedLogLike;
	/** loglikelihood after this step
	 * ({@link #accepted} ? {@link #proposedLogLike} : {@link #oldLogLike}) */
	public double newLogLike;
	/** heat of the chain */
	public double heat;
	/** log of (backproposal probability / proposal probability) */
	public double bpp;
	
	/** <code>true</code> if this step was accepted */
	public boolean accepted;
	
	// TODO add step type
	
	// TODO add statistics (acceptance rates etc)
}
