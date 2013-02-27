package statalign.base.mcmc;

import statalign.base.mcmc.PriorDistribution;
import statalign.model.ext.plugins.structalign.Funcs;
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
		return Funcs.logGammaDensity(x, shape, rate);
	}
}
