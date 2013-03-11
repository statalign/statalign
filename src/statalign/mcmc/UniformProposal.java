package statalign.mcmc;

import statalign.base.Utils;
public class UniformProposal extends ProposalDistribution<Double> {

	private double centre;
	private double width;
	
	public double logDensity(Double x) {
		if ((x > (centre-width/2)) && (x < (centre+width/2))) {
			return 0.0;
		}
		else {
			return Double.NEGATIVE_INFINITY;
		}
	}
	public Double sample() {
		return (Utils.generator.nextDouble() - 0.5) * width + centre;
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		centre = currentParam;
		width = proposalWidthControlVariable;
	}
	
}
