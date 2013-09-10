package statalign.mcmc;

import org.apache.commons.math3.util.FastMath;
import statalign.utils.NormalDistribution;

public class LogNormalPrior implements PriorDistribution<Double> {
	
	private NormalDistribution n;
	private double mean;
	private double sd;
	
	public LogNormalPrior(double a, double b) {
		mean = a;
		sd = b;
		n = new NormalDistribution(mean,sd);
	}
	public double logDensity(Double x) {
		return Math.log(n.density(FastMath.log(x)));
	}
	public double logDensityUnnormalised(Double x) {
		return -Math.pow((FastMath.log(x)-mean)/sd, 2)/2 - FastMath.log(sd);
	}
}
