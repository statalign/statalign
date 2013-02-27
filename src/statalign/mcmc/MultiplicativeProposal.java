package statalign.mcmc;

import statalign.utils.NormalDistribution;
public class MultiplicativeProposal extends ProposalDistribution<Double> {

	private NormalDistribution n;
	
	public MultiplicativeProposal() {
		n = new NormalDistribution(0,1);
	}
	public double logDensity(Double x) {
		return x * Math.log(n.density(x));
	}
	public Double sample() {
		return Math.exp(n.sample());
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		n = new NormalDistribution(Math.log(currentParam),proposalWidthControlVariable);
	}
	
}
