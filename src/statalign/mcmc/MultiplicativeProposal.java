package statalign.mcmc;

import statalign.base.Utils;
import statalign.utils.NormalDistribution;
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
		if (Utils.DEBUG) System.out.println("logDensity(x) of "+x);
		//return Math.log(n.density(x)) - Math.log(x);
		return n.logDensity(x) - Math.log(x);
	}
	public Double sample() {
		return Math.exp(n.sample());
	}
	public void updateProposal(double proposalWidthControlVariable, 
			Double currentParam) {
		if (Utils.DEBUG) System.out.println("New currentParam = "+currentParam+", var = "+proposalWidthControlVariable);
		n = new NormalDistribution(Math.log(currentParam),proposalWidthControlVariable);
	}
	
}
