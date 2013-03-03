package statalign.mcmc;

import statalign.base.Utils;
import statalign.mcmc.PriorDistribution;
import statalign.utils.GammaDistribution;


public class GammaPrior implements PriorDistribution<Double> {
	
	private GammaDistribution g;
	private double shape;
	private double rate;
	
	public GammaPrior(double a, double b) {
		shape = a;
		rate = b;
		g = new GammaDistribution(shape,1.0/rate);
	}
	public double logDensity(Double x) {
		return Math.log(g.density(x));
	}
	public double logDensityUnnormalised(Double x) {
		return Utils.logGammaDensity(x, shape, rate);
	}
}
