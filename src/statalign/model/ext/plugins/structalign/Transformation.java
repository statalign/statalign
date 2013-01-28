package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import statalign.model.ext.plugins.structalign.Funcs;
import statalign.utils.NormalDistribution;

public class Transformation{
	//RealVector axis;
	//double rot;
	RealVector xlat;
	//RealMatrix rotMatrix;
	Rotation rotMatrix;
	// sum of squares for coordinate sub-matrix optimal rotation
	double ss;
	
	Transformation(){}
	
	Transformation(Vector3D axis, double angle){
		rotMatrix = new Rotation(axis, angle);
	}
	
	Transformation(double[] axis, double angle, double[] xlat){
		Vector3D ax = new Vector3D(axis[0], axis[1], axis[2]);
		rotMatrix = new Rotation(ax, angle);
		this.xlat = new ArrayRealVector(xlat);
	}
	
	Transformation(Vector3D axis, double angle, double[] xlat){
		rotMatrix = new Rotation(axis, angle);
		this.xlat = new ArrayRealVector(xlat);
	}
	
	public RealMatrix getRealRotation(){
		return new Array2DRowRealMatrix(rotMatrix.getMatrix());
	}
	
	/**
	 * 
	 * @return same Transformation object with axis and angle filled
	 */		
	/* public Transformation fillAxisAngle(){
		EigenDecomposition eigen = new EigenDecomposition(rotMatrix);
		
		double real = 0;
		int i = -1;
		// TODO seems like there should be a better way to do this
		while(real < .99999999){
			i++;
			real = eigen.getRealEigenvalue(i);
		}
		axis = eigen.getEigenvector(i);
		
		RealVector v;
		if(axis.getEntry(2) < .98)
			v = Funcs.crossProduct(axis, new ArrayRealVector(new double[] {0, 0, 1}));
		else
			v = Funcs.crossProduct(axis, new ArrayRealVector(new double[] {1, 0, 0}));
		
		RealVector rotv = rotMatrix.preMultiply(v);
		
		rot = Math.acos(v.dotProduct(rotv) / (v.getNorm() * rotv.getNorm()) );
		// check axis
		axis = axis.mapMultiply(Math.signum(axis.dotProduct(Funcs.crossProduct(v, rotv))));
		return this;
	} */
	
	/** 
	 * Calculates the rotation matrix from an axis and angle
	 * with axis x and angle r, the rotation matrix is
	 * cos(r)I + sin(r)cross(x) + (1-cos(r))outer(x)
	 * where cross(x) is the crossproduct matrix of x and outer(x) is the outer (tensor) product
	 * 
	 * @return Transformation object with rotation matrix calculated from axis and angle
	 */
	/* public Transformation fillRotationMatrix(){
		RealMatrix outer = axis.outerProduct(axis);
		// use transpose of cross product matrix (as given on wikipedia) because
		// we post-mulitply by rotation matrix (other components of final step are symmetric)
		RealMatrix crossTranspose = new Array2DRowRealMatrix(new double[][] { 
				{0, -axis.getEntry(2), axis.getEntry(1)},
				{axis.getEntry(2), 0, -axis.getEntry(0)},
				{-axis.getEntry(1), axis.getEntry(0), 0} 
		});
		RealMatrix Icos = new Array2DRowRealMatrix(new double[3][3]);
		for(int i = 0; i < 3; i++)
			Icos.setEntry(i, i, Math.cos(rot));
		crossTranspose = crossTranspose.scalarMultiply(Math.sin(rot));
		outer = outer.scalarMultiply(1 - Math.cos(rot));
		rotMatrix = Icos.add(crossTranspose).add(outer);
		return this;
	} */
	
	public double density(Transformation candidate, double kvM, double kvMF, double sd){
		double density = vonMises.density(kvM, rotMatrix.getAngle(), candidate.rotMatrix.getAngle());
		density *= vonMisesFisher.density(kvMF, rotMatrix.getAxis(), candidate.rotMatrix.getAxis());
		for(int i = 0; i < 3; i++)
			density *= new NormalDistribution(xlat.getEntry(i), sd).density(candidate.xlat.getEntry(i));
		return density;
	}
	
	public void lsTrans(RealMatrix A, RealMatrix B){
		
		int n = A.getRowDimension();
		
		RealMatrix C = new Array2DRowRealMatrix(new double[n][n]);
		for(int i = 0; i < n; i++)
			C.setEntry(i, i, 1);
		C = C.scalarAdd(-1/ (double) n);
		
		SingularValueDecomposition svd = new SingularValueDecomposition(B.transpose().multiply(C).multiply(A));

		// U V^T from the SVD yields the optimal rotation/reflection matrix
		// for rotations only, use U S V^T, where S is the identity with
		// S_{n,n} = det(U V^T)
		double det = new LUDecomposition(svd.getU().multiply(svd.getVT())).getDeterminant();
		RealMatrix S = new Array2DRowRealMatrix(new double[3][3]);
		S.setEntry(0, 0, 1);
		S.setEntry(1, 1, 1);
		S.setEntry(2, 2, det);

		// rotation = U S V^T  
		rotMatrix = new Rotation( svd.getU().multiply(S).multiply(svd.getVT()).getData(), 1e-20 );
		
		RealMatrix temp = getRealRotation();
		
		// translation = mean(a) - mean(b) R
		xlat = Funcs.meanVector(A).subtract(temp.preMultiply(Funcs.meanVector(B))) ;
		
		// with Rotation class, use transpose (revert method)
		rotMatrix = rotMatrix.revert();		
	}
	
	// </Transformation>
}