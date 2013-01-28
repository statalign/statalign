package statalign.model.ext.plugins.structalign;

import statalign.base.Utils;
import cern.jet.math.Bessel;

public class vonMises{
	/**
	 * 
	 * @param kappa - concentration parameter
	 * @param mean - mean angle
	 * @param angle - angle to be evaluated
	 * @return log density of distribution at @param angle
	 */
	static double logDensity(double kappa, double mean, double angle){
		return kappa * Math.cos(angle - mean) - Math.log(2 * Math.PI * Bessel.i0(kappa));
	}
	
	static double density(double kappa, double mean, double angle){
		return Math.exp(kappa * Math.cos(angle - mean)) / (2 * Math.PI * Bessel.i0(kappa));
	}
	
	/**
	 * 
	 * @param kappa - concentration parameter
	 * @param mean - mean angle
	 * @return - angle simulated from vonMises distribution with given parameters
	 */
	static double simulate(double kappa, double mean){
		double vm = 0, U1, U2, U3, a, b, c, f, r, z;
		boolean cont;
		a = 1 + Math.pow(1 + 4 * Math.pow(kappa, 2), 0.5);
		b = (a - Math.pow(2 * a, 0.5)) / (2 * kappa);
		r = (1 + Math.pow(b, 2)) / (2 * b);
		cont = true;
		while (cont) {
			U1 = Utils.generator.nextDouble();
			z = Math.cos(Math.PI * U1);
			f = (1 + r * z) / (r + z);
			c = kappa * (r - f);
			U2 = Utils.generator.nextDouble();
			if (c * (2 - c) - U2 > 0) {
				U3 = Utils.generator.nextDouble();
				vm = Math.signum(U3 - 0.5) * Math.acos(f) + mean;
				vm = vm % (2 * Math.PI);
				cont = false;
			}
			else {
				if (Math.log(c/U2) + 1 - c >= 0) {
					U3 = Utils.generator.nextDouble();
					vm = Math.signum(U3 - 0.5) * Math.acos(f) + mean;
					vm = vm % (2 * Math.PI);
					cont = false;
				}
			}
		}
		return vm;
	}
	// </vonMises>
}
