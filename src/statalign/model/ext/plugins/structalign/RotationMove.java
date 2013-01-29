package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import statalign.base.Utils;
import statalign.model.ext.plugins.StructAlign;
import statalign.model.ext.plugins.structalign.vonMises;

public class RotationMove extends RotationOrTranslationMove {
	
	public RotationMove (StructAlign s, String n) {
		owner = s;
		name = n;
		proposalWidthControlVariable = 1.0/owner.angleP;
	}

	public double proposal(Object externalState) {
		double smallAngle = vonMises.simulate(1.0/proposalWidthControlVariable, 0);
		
		RealVector randomAxis = new ArrayRealVector(3);
		for(int i = 0; i < 3; i++)
			randomAxis.setEntry(i, Utils.generator.nextGaussian());
		randomAxis.unitize();
		
		Rotation Q = new Rotation(new Vector3D(randomAxis.toArray()), smallAngle);
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			Rotation R = new Rotation(new Vector3D(owner.axes[j]), owner.angles[j]);
		
			R = Q.applyTo(R);
			
			owner.axes[j] = R.getAxis().toArray();
			owner.angles[j] = R.getAngle();
		}
					
		return 0;
		// logProposalRatio is 0 because prior is uniform and proposal is symmetric
	}	
}
