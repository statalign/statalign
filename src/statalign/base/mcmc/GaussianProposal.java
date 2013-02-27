package statalign.base.mcmc;

import statalign.utils.NormalDistribution;
public class GaussianProposal extends ProposalDistribution<Double> {

	private NormalDistribution n;
	
	public GaussianProposal() {
		n = new NormalDistribution(0,1);
	}
	public double logDensity(Double x) {
		return Math.log(n.density(x));
	}
	public Double sample() {
		return n.sample();
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		n = new NormalDistribution(currentParam,proposalWidthControlVariable);
	}
	
}
