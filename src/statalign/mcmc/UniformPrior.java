package statalign.mcmc;


public class UniformPrior implements PriorDistribution<Double> {

	double min = 0.0;
	double max = Double.POSITIVE_INFINITY;
	
	public UniformPrior() { }
	public UniformPrior(double a, double b) { 
		if (a >= b) {
			throw new IllegalArgumentException("Lower bound, "+a+", in UniformPrior is not below the upper bound, "+b+".");
		}
		min = a; 
		max = b;
	}
	public double logDensity(Double x) {
		if (max == Double.POSITIVE_INFINITY) {
			throw new RuntimeException("Cannot normalise density for UniformPrior.");
		}
		else {
			return (logDensityUnnormalised(x) - Math.log(max-min));
		}
	}
	public double logDensityUnnormalised(Double x) {
		if (x > min && x < max) {
			return 0;
		}
		else {
			return Double.NEGATIVE_INFINITY;
		}
	}
}
