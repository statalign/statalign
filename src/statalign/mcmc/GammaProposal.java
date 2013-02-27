package statalign.mcmc;

import statalign.utils.GammaDistribution;
public class GammaProposal extends ProposalDistribution<Double> {

	private GammaDistribution g;
	private double proposalShape;
	private double proposalRate;
	
	public GammaProposal(double a, double b) {
		proposalShape = a;
		proposalRate = b;
		g = new GammaDistribution(a,1.0/b);
	}
	public double logDensity(Double x) {
		return Math.log(g.density(x));
	}
	public Double sample() {
		return g.sample();
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		double scale = currentParam < 1e-6 ? 1e-6 : currentParam;
		double conc = 1.0/proposalWidthControlVariable;
		g = new GammaDistribution(conc + proposalShape, 
		scale / (conc + proposalRate));
	}
	
}
