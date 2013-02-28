package statalign.mcmc;

import statalign.utils.NormalDistribution;
public class LogisticProposal extends ProposalDistribution<Double> {

	private NormalDistribution n;
	
	public LogisticProposal() {
		n = new NormalDistribution(0,1);
	}
	/**
	 * Computes the logDensity of proposing a value x, adding
	 * on the term from the Jacobian of the transformation from
	 * x to logit(x).
	 */
	public double logDensity(Double x) {
		return Math.log(n.density(x)) - Math.log(x*(1-x));
	}
	@Override
	public Double sample() { 
		double newSample = Math.exp(n.sample());
		return (newSample / (1 + newSample));
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		n = new NormalDistribution(Math.log(currentParam) - 
									Math.log(1-currentParam),
				                      proposalWidthControlVariable);
	}
	
}
