package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import statalign.base.Utils;
import statalign.model.ext.plugins.StructAlign;
import statalign.model.ext.plugins.structalign.Transformation;
import statalign.model.ext.plugins.structalign.Funcs;


public class RotationProposal{
	
	StructAlign owner;
	/* libraries of axes, angles, and translation vectors for each protein */
	Transformation[][] libraries;
	int window = 25;
	int fixed = 0;
	double percent = 1;  // (value between 0 and 100)
	
	/* concentration parameters of proposal distribution */
	double kvMF = 100;
	double kvM = 100;
	public double sd = .1;
	
	// TODO All of the values above this point should probably be chosen elsewhere
	
	public RotationProposal(StructAlign s){
		owner = s;
		// calculate all rotations relative to fixed protein
		libraries = new Transformation[owner.coords.length][];
		for(int i = 0; i < owner.coords.length; i++)
			if(i != fixed)
				libraries[i] = selectBest(calculateAllOptimal(fixed, i, window), percent);
	}
	
	public Transformation[] calculateAllOptimal(int a, int b, int window){
		RealMatrix A = new Array2DRowRealMatrix(owner.coords[a]);
		RealMatrix B = new Array2DRowRealMatrix(owner.coords[b]);
		int nA = A.getColumn(0).length - window + 1;
		int nB = B.getColumn(0).length - window + 1;
		Transformation[] allOptimal = new Transformation[nA * nB];
		for(int i = 0; i < nA * nB; i++)
			allOptimal[i] = new Transformation();
		
		// create centering matrix 
		RealMatrix C = new Array2DRowRealMatrix(new double[window][window]);
		for(int i = 0; i < window; i++)
			C.setEntry(i, i, 1);
		C = C.scalarAdd(-1/ (double) window);
		
		RealMatrix subA;
		RealMatrix subB;
		
		int k = 0;
		for(int i = 0; i < nA; i++){
			for(int j = 0; j < nB; j++){
				subA = A.getSubMatrix(i, i + window - 1, 0, 2);
				subB = B.getSubMatrix(j, j + window - 1, 0, 2);
				
				SingularValueDecomposition svd = new SingularValueDecomposition(subB.transpose().multiply(C).multiply(subA));
				
				// U V^T from the SVD yields the optimal rotation/reflection matrix
				// for rotations only, use U S V^T, where S is the identity with
				// S_{n,n} = det(U V^T)
				double det = new LUDecomposition(svd.getU().multiply(svd.getVT())).getDeterminant();
				RealMatrix S = new Array2DRowRealMatrix(new double[3][3]);
				S.setEntry(0, 0, 1);
				S.setEntry(1, 1, 1);
				S.setEntry(2, 2, det);
				
				// rotation = U S V^T
				RealMatrix R = svd.getU().multiply(S).multiply(svd.getVT());
				
				// translation = mean(a) - mean(b) R
				RealVector xlat = Funcs.meanVector(subA).subtract(R.preMultiply(Funcs.meanVector(subB)));
				
				// calculate the resulting sum of squares
				double ss = 0;
				
				// a - (R b + 1 xlat)
				RealVector one = new ArrayRealVector(new double[window]).mapAdd(1);
				
				double[][] diff = subA.subtract( subB.multiply(R).add(one.outerProduct(xlat))).getData();
				for(int m = 0; m < diff.length; m++)
					for(int n = 0; n < diff[0].length; n++)
						ss += Math.pow(diff[m][n], 2);
				allOptimal[k].ss = ss;
				allOptimal[k].rotMatrix = new Rotation(R.getData(), 1e-20).revert(); // Rotation class applies rotations
				allOptimal[k].xlat = xlat;											 // differently, so use inverse
				k++;
			}
		}	
		return allOptimal;
	}
	
	/**
	 * @param
	 * @param library - a set of optimal Transformations between protein fragments
	 * @return - the @percent of optimal Transformations with the lowest RMSD
	 */
	
	public Transformation[] selectBest(Transformation[] library, double percent){
		double[] ss = new double[library.length];
		for(int i = 0; i < ss.length; i++)
			ss[i] = library[i].ss;
		double cutoff = new Percentile().evaluate(ss, percent);
		int n = 0;
		for(int i = 0; i < library.length; i++)
			n += library[i].ss < cutoff ? 1 : 0;
		Transformation[] best = new Transformation[n];
		int j = 0;
		for(int i = 0; i < library.length; i++)
			if(library[i].ss < cutoff)
				best[j++] = library[i];
		return best;
	}
	

	/** propose from a library mixture distribution */
	public Transformation propose(int index){
		Transformation[] library = libraries[index];
		int n = library.length;
		int k = Utils.generator.nextInt(n);
		
		Vector3D axis = vonMisesFisher.simulate(kvMF, library[k].rotMatrix.getAxis());
		double angle = vonMises.simulate(kvM, library[k].rotMatrix.getAngle());
		double[] xlat = new double[3];
		for(int i = 0; i < 3; i++)
			xlat[i] = Utils.generator.nextGaussian() * sd + library[k].xlat.getEntry(i);			
		return new Transformation(axis, angle, xlat);
	}
	
	/**
	 * 
	 * @param index - index indicating the set of Transformations on which mixture components are centered
	 * @param candidate - a candidate Transformation whose density is to be evaluated
	 * @return the log density of @candidate according to the library mixture distribution
	 */
	public double libraryLogDensity(int index, Transformation candidate){
		Transformation[] library = libraries[index];
		double density = 0;
		for(int i = 0; i < library.length; i++)
			density += library[i].density(candidate, kvM, kvMF, sd);
		density /= library.length;
		return Math.log(density);
	}						
	// </RotationProposal>
}