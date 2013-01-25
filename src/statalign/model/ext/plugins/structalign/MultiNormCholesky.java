package statalign.model.ext.plugins.structalign;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;

/** Adapted from org.apache.commons.math3.distribution.MultivariateNormalDistribution
 * 
 * @author Challis
 * 
 */

public class MultiNormCholesky{
	/** Dimension. */
	private final int dim;
	/** Vector of means. */
	private final double[] means;
	/** Covariance matrix. */
	private final RealMatrix covarianceMatrix;
	/** The matrix inverse of the covariance matrix. */
	private final RealMatrix covarianceMatrixInverse;
	/** The determinant of the covariance matrix. */
	private double covarianceMatrixDeterminant;
	/** Matrix used in computation of samples. */

	/**
	 * Creates a multivariate normal distribution with the given mean vector and
	 * covariance matrix.
	 * <br/>
	 * The number of dimensions is equal to the length of the mean vector
	 * and to the number of rows and columns of the covariance matrix.
	 * It is frequently written as "p" in formulae.
	 *
	 * @param means Vector of means.
	 * @param covariances Covariance matrix.
	 * @throws DimensionMismatchException if the arrays length are
	 * inconsistent.
	 * @throws SingularMatrixException if the eigenvalue decomposition cannot
	 * be performed on the provided covariance matrix.
	 */
	public MultiNormCholesky(final double[] means,
			final double[][] covariances)
					throws SingularMatrixException,
					DimensionMismatchException,
					NonPositiveDefiniteMatrixException {

		dim = means.length;

		if (covariances.length != dim) {
			throw new DimensionMismatchException(covariances.length, dim);
		}

		for (int i = 0; i < dim; i++) {
			if (dim != covariances[i].length) {
				throw new DimensionMismatchException(covariances[i].length, dim);
			}
		}

		this.means = MathArrays.copyOf(means);

		covarianceMatrix = new Array2DRowRealMatrix(covariances);

		// Covariance matrix eigen decomposition.
		final CholeskyDecomposition covMatDec;
		try {
			covMatDec = new CholeskyDecomposition(covarianceMatrix);
		}
		catch (NonPositiveDefiniteMatrixException e) {
			System.out.println(e);
			System.out.println("");
			System.out.println("covariances = ");
			for (int i=0; i<dim; i++) {
				for (int j=0; j<dim; j++) {
					System.out.print(" "+covariances[i][j]);
				}
				System.out.println("");
			}	
			throw new RuntimeException(e);
		}

		// Compute and store the inverse.
		covarianceMatrixInverse = covMatDec.getSolver().getInverse();
		// Compute and store the determinant.
		covarianceMatrixDeterminant = 0;
		for(int i = 0; i < dim; i++)
			covarianceMatrixDeterminant += 2 * Math.log(covMatDec.getL().getEntry(i, i));			
	}

	/**
	 * Gets the mean vector.
	 *
	 * @return the mean vector.
	 */
	public double[] getMeans() {
		return means;
	}

	/**
	 * Gets the covariance matrix.
	 *
	 * @return the covariance matrix.
	 */
	public RealMatrix getCovariances() {
		return covarianceMatrix;
	}

	/** {@inheritDoc} */
	public double logDensity(final double[] vals) throws DimensionMismatchException {
		if (vals.length != dim) {
			throw new DimensionMismatchException(vals.length, dim);
		}

		double x = (double)-dim / 2 * FastMath.log(2 * FastMath.PI) + 
				// -0.5 * FastMath.log(covarianceMatrixDeterminant) +
				-0.5 * covarianceMatrixDeterminant +
				getExponentTerm(vals);

		/*System.out.println(dim);
			System.out.println((double)-dim / 2 * FastMath.log(2 * FastMath.PI));
			System.out.println(-0.5 * covarianceMatrixDeterminant);
			System.out.println(getExponentTerm(vals));*/

		return x;
	}

	/**
	 * Computes the term used in the exponent (see definition of the distribution).
	 *
	 * @param values Values at which to compute density.
	 * @return the multiplication factor of density calculations.
	 */
	private double getExponentTerm(final double[] values) {
		final double[] centered = new double[values.length];
		for (int i = 0; i < centered.length; i++) {
			centered[i] = values[i] - getMeans()[i];
		}
		final double[] preMultiplied = covarianceMatrixInverse.preMultiply(centered);
		double sum = 0;
		for (int i = 0; i < preMultiplied.length; i++) {
			sum += preMultiplied[i] * centered[i];
		}
		return -0.5 * sum;
	}
}
