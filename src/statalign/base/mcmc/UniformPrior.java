package statalign.base.mcmc;


public class UniformPrior implements PriorDistribution<Double> {

	public double logDensity(Double x) {
		throw new RuntimeException("Cannot normalise density for UniformPrior.");
	}
	public double logDensityUnnormalised(Double x) {
		return 0;
	}
}
