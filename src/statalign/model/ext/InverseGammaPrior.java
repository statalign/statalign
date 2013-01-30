package statalign.model.ext;

import statalign.model.ext.PriorDistribution;
import statalign.model.ext.plugins.structalign.Funcs;
import statalign.utils.GammaDistribution;

public class InverseGammaPrior implements PriorDistribution<Double> {
	
	private GammaDistribution g;
	private double shape;
	private double rate;
	
	public InverseGammaPrior(double a, double b) {
		shape = a;
		rate = b;
		g = new GammaDistribution(shape,1.0/rate);
	}
	public double logDensity(Double x) {
		return Math.log(g.density(1.0/x));
	}
	public double logDensityUnnormalised(Double x) {
		return Funcs.logGammaDensity(1.0/x, shape, rate);
	}
}
