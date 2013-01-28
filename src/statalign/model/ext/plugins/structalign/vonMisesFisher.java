package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import statalign.base.Utils;
import statalign.utils.BetaDistribution;

public class vonMisesFisher{
	/**
	 * 
	 * @param kappa - concentration parameter
	 * @param mean - mean direction
	 * @param v - vector to be evaluated
	 * @return log density of distribution at @param v
	 */
	static double logDensity(double kappa, Vector3D mean, Vector3D v){
		return Math.log(kappa) - Math.log(4 * Math.PI * Math.sinh(kappa)) + kappa * mean.dotProduct(v);
	}	
	
	static double density(double kappa, Vector3D mean, Vector3D v){
		return kappa / (4 * Math.PI * Math.sinh(kappa)) * Math.exp(kappa * mean.dotProduct(v));
	}	
	
	/**
	 * 
	 * @param kappa - concentration parameter
	 * @param v - mean vector
	 * @return - vector simulated from von Mises-Fisher distribution with given parameters
	 */
	// TODO replace with exact simulation method for S2 from Kent et al 2012
	static Vector3D simulate(double kappa, Vector3D mean){
		double b, x0, c, z, u, w;
		int dim = 3;
		RealVector unitSphere = new ArrayRealVector(new double[dim-1]);
		RealVector sample;
		do {
			b = (-2*kappa + Math.pow(Math.pow(4*Math.pow(kappa,2) + (dim-1),2),.5) ) / (dim-1);
			x0 = (1-b)/(1+b);
			c = kappa*x0 + (dim-1)*Math.log(1-Math.pow(x0, 2));
			z = new BetaDistribution((dim-1)/2,(dim-1)/2).sample();
			u = Utils.generator.nextDouble();
			w = (1 - (1+b)*z)/(1 - (1-b)*z);
		} while(kappa*w + (dim-1) * Math.log(1-x0*w) - c < Math.log(u));

		for(int i = 0; i < dim - 1; i++)
			unitSphere.setEntry(i, Utils.generator.nextGaussian());
		unitSphere = unitSphere.unitVector();
		sample = unitSphere.mapMultiply(Math.pow(1-Math.pow(w,2), .5));
		Vector3D sample3D = new Vector3D(sample.getEntry(0), sample.getEntry(1), w);
				  
		// calculate rotation to rotate (0,0,1) to mean
		Rotation Q = new Rotation(new Vector3D(0, 0, 1), mean);
		
		// Output simulated vector rotated to center around mean
		sample3D = Q.applyTo(sample3D);
		
		return sample3D;
	}
		// </vonMisesFisher>
}
