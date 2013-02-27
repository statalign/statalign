package statalign.mcmc;

import statalign.utils.NormalDistribution;
public class LogisticProposal extends ProposalDistribution<Double> {

	private NormalDistribution n;
	
	public LogisticProposal() {
		n = new NormalDistribution(0,1);
	}
	public double logDensity(Double x) {
		return x * (1-x) * Math.log(n.density(x));
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
