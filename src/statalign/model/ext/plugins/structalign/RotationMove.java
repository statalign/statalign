package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import statalign.base.Utils;
import statalign.model.ext.plugins.StructAlign;
import statalign.model.ext.plugins.structalign.vonMises;

public class RotationMove extends RotationOrTranslationMove {

	StructAlign owner;
	
	public RotationMove (StructAlign s) {
		owner = s;
		omit = 1;
	}

	public double proposal() {
		double smallAngle = vonMises.simulate(owner.angleP, 0);
		
		RealVector randomAxis = new ArrayRealVector(3);
		for(int i = 0; i < 3; i++)
			randomAxis.setEntry(i, Utils.generator.nextGaussian());
		randomAxis.unitize();
		
		Rotation Q = new Rotation(new Vector3D(randomAxis.toArray()), smallAngle);
		Rotation R = new Rotation(new Vector3D(owner.axes[ind]), owner.angles[ind]);
		
		R = Q.applyTo(R);
	
		owner.axes[ind] = R.getAxis().toArray();
		owner.angles[ind] = R.getAngle();
				
		return 0;
		// logProposalRatio is 0 because prior is uniform and proposal is symmetric
	}
}
