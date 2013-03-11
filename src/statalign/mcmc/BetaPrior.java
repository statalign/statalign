package statalign.mcmc;

import statalign.base.Utils;
import statalign.mcmc.PriorDistribution;
import statalign.utils.BetaDistribution;


public class BetaPrior implements PriorDistribution<Double> {
	
	private BetaDistribution b;
	private double alpha;
	private double beta;
	
	public BetaPrior(double _alpha, double _beta) {
		alpha = _alpha;
		beta = _beta;
		b = new BetaDistribution(alpha,beta);
	}
	public double logDensity(Double x) {
		return Math.log(b.density(x));
	}
	public double logDensityUnnormalised(Double x) {
		return Utils.logBetaDensity(x, alpha, beta);
	}
}
