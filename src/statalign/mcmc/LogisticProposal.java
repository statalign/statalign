package statalign.mcmc;

import statalign.base.Utils;
import statalign.utils.NormalDistribution;
/**
 * Proposal distribution for variables whose support is [0,1].
 * In many cases this type of proposal results in higher effective
 * sample size and better tail coverage than, for example, a simple
 * Gaussian random walk. 
 * 
 * @author herman 
 */
public class LogisticProposal extends ProposalDistribution<Double> {

	private NormalDistribution n;
	
	private double logit(double x) {
		return Math.log(x) - Math.log(1-x);
	}
	
	public LogisticProposal() {
		n = new NormalDistribution(0,1);
	}
	/**
	 * Computes the logDensity of proposing a value x, adding
	 * on the term from the Jacobian of the transformation from
	 * x to logit(x).
	 */
	public double logDensity(Double x) {
		return n.logDensity(logit(x)) - Math.log(x*(1-x));
	}
	@Override
	public Double sample() { 
		double newSample = Math.exp(n.sample());
		return (newSample / (1 + newSample));
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		n = new NormalDistribution(logit(currentParam),
				                      proposalWidthControlVariable);
	}
	
}
