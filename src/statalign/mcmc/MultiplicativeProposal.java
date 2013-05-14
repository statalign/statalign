package statalign.mcmc;

import statalign.base.Utils;
import statalign.utils.NormalDistribution;
/**
 * Proposal that multiplies the current parameter by a log-normally
 * distributed variable, which is equivalent to a Gaussian random
 * walk on the logarithm of the parameter. 
 * 
 * @author herman
 */
public class MultiplicativeProposal extends ProposalDistribution<Double> {

	private NormalDistribution n;
	
	public MultiplicativeProposal() {
		n = new NormalDistribution(0,1);
	}
	/**
	 * Computes the logDensity of proposing a value x, adding
	 * on the term from the Jacobian of the transformation from
	 * x to log x.
	 */
	public double logDensity(Double x) {
		return n.logDensity(Math.log(x)) - Math.log(x);
	}
	public Double sample() {
		return Math.exp(n.sample());
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		n = new NormalDistribution(Math.log(currentParam),proposalWidthControlVariable);
	}		
	
}
