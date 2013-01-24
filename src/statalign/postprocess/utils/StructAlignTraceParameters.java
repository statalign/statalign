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
	 * The parameter values
	 */
	public double[] sigma2;
	public boolean[] sigma2Proposed;
	public double sigma2Hier;
	public boolean sigma2HProposed;
	public double nu;
	public boolean nuProposed;
	public double tau;
	public boolean tauProposed;
	public double epsilon;
	public boolean epsilonProposed;
	
	public boolean globalSigma = true;
	/**
	 * Whether the sample comes from the burn-in phase or not.
	 */
	public boolean burnin;
	
	/**
	 * This constructor initialises the container.
	 * @param burnin True if the sample comes from the burn-in phase.
	 */
	public StructAlignTraceParameters(boolean burnin){
		this.burnin = burnin;
//		
	}
}
