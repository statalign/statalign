package statalign.model.ext;

import statalign.model.ext.PriorDistribution;
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
	public void updateDistribution(double a, double b) {
		g = new GammaDistribution(shape + a,1.0/(b + rate));
	}
	public double logDensity(Double x) {
		return Math.log(g.density(x));
	}
}
