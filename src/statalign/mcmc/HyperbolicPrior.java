package statalign.mcmc;

import statalign.mcmc.PriorDistribution;
import statalign.model.ext.plugins.structalign.Funcs;
import statalign.utils.GammaDistribution;


public class HyperbolicPrior implements PriorDistribution<Double> {

	public double logDensity(Double x) {
		if (false) {
				return 0;
		}
		else {
			throw new RuntimeException("Cannot normalise density for HyperbolicPrior.");
		}		
	}
	public double logDensityUnnormalised(Double x) {
		return 1.0/Math.sqrt(x);
	}
}
