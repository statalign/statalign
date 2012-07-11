package statalign.model.subst.plugins;

import statalign.model.subst.plugins.NucleotideModel;
import statalign.base.*;
import statalign.Jama.*; 

/**
 * This class implements various utility methods for nucleotide models.
 * 
 * @author lyngsoe
 */
public class NucleotideModelUtils {

	/**
	 * Class method for sampling new rate parameter for one of the
	 * rates in r, where the maximum allowed change is dmax. Return
	 * value is logarithm of Metropolis-Hastings ratio.
	 */
	public static double sampleRates(double r[], double dmax){
		int k = Utils.generator.nextInt(r.length);
		double l = Math.max(0, r[k] - dmax);
		double span = r[k] + dmax - l;
		do{ /* Rates of zero are disallowed */
			r[k] = l + Utils.generator.nextDouble() * span;
		} while (r[k] <= 0);
		return Math.log((r[k] + dmax - Math.max(0, r[k] - dmax)) / span);
	}

	/**
	 * Class method for sampling new equilibrium frequency parameter
	 * for one of the frequencies in f, where the maximum allowed
	 * change is dmax. Return value is logarithm of
	 * Metropolis-Hastings ratio.
	 */
	public static double sampleFrequencies(double f[], double dmax){
		int k = Utils.generator.nextInt(f.length);
		double l = Math.max(0, f[k] - dmax);
		double r = Math.min(1, f[k] + dmax);
		double x;
		boolean resample = true;
		do{ /* Frequencies are not allowed to be zero or one */
			do{
				x = l + Utils.generator.nextDouble() * (r - l);
			} while ((x <= 0) || (x >= 1));
			resample = false;
			for (int i = 0; i < f.length; i++)
				if (i != k){
					f[i] *= (1.0 - x) / (1.0 - f[k]);
					if ((f[i] <= 0) || (f[i] >= 1))
						resample = true;
				}
			f[k] = x;
		} while (resample);
		return Math.log((Math.min(1.0, f[k] + dmax) - Math.max(0.0, f[k] - dmax)) / (r - l));
	}

	/**
	 * Class method for numerical diagonalisation of rate matrix Q for
	 * nucleotide model M. It is assumed that M.e holds the
	 * equilibrium frequencies for the model and that the
	 * diagonalising matrices should be stored in M.v, M.d and M.w,
	 * and that it is safe to modify Q.
	 */
	public static void numericalDiagonalisation(NucleotideModel M, double Q[][]){
		/* Symmetricalise Q */
		for (int i = 0; i < Q.length; i++)
			for (int j = 0; j < i; j++)
				Q[j][i] = Q[i][j] = Math.sqrt((double)M.e[i] / M.e[j]) * Q[i][j];

		/* Find matrix of eigenvectors */
		Matrix R = new Matrix(Q);
		EigenvalueDecomposition E = new EigenvalueDecomposition(R);
		Matrix V = E.getV();
		double d[] = E.getRealEigenvalues();

		/* Copy information to M */
		for (int i = 0; i < Q.length; i++){
			M.d[i] = d[i];
			for (int j = 0; j < Q.length; j++){
				M.v[i][j] = Math.sqrt(1.0 / M.e[i]) * V.get(i, j);
				M.w[j][i] = Math.sqrt(M.e[i]) * V.get(i, j);
			}
		}
	}
}

