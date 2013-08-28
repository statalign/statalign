package statalign.mcmc;

import statalign.base.Utils;
import statalign.mcmc.PriorDistribution;
import statalign.utils.GammaDistribution;

public class LinearPrior implements PriorDistribution<Double> {
	
	private double gradient;
	
	public LinearPrior() { gradient = 1; }		
	public LinearPrior(double a) {
		gradient = a;
	}
	public double logDensity(Double x) {
		throw new RuntimeException("Cannot normalise improper linear density.");
	}
	public double logDensityUnnormalised(Double x) {
		return Math.log(x) + Math.log(gradient);
	}
}
