package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.linear.ArrayRealVector;

import statalign.model.ext.plugins.StructAlign;
import statalign.model.ext.plugins.structalign.Transformation;

public class LibraryMove extends RotationOrTranslationMove {
	
	public LibraryMove (StructAlign s, String n) {
		owner = s;
		name = n;
	}

	public double proposal(Object externalState) {
		
		owner.rotProp.sd = proposalWidthControlVariable;
		Transformation oldSub = new Transformation(owner.axes[index], owner.angles[index], owner.xlats[index]);
		// transformation should be relative to reference protein
		oldSub.xlat = oldSub.xlat.subtract(new ArrayRealVector(owner.xlats[0]));
		Transformation libProp = owner.rotProp.propose(index);
		owner.axes[index] = libProp.rotMatrix.getAxis().toArray();
		owner.angles[index] = libProp.rotMatrix.getAngle();
		owner.xlats[index] = libProp.xlat.toArray();

		// library density 
		double logProposalRatio = owner.rotProp.libraryLogDensity(index, oldSub) - 
				owner.rotProp.libraryLogDensity(index, libProp);

		// proposed translation is relative to reference protein
		for(int i = 0; i < 3; i++)
			owner.xlats[index][i] += owner.xlats[0][i];			

		// calculate 'difference' between proposed and current transformations
		double[] diffxlat = new double[3];
		for(int i = 0; i < 3; i++)
			diffxlat[i] = owner.xlats[index][i] - oldxlats[index][i];
		Rotation diffRotMat = libProp.rotMatrix.applyTo(oldSub.rotMatrix.revert());

		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			if(j != index){
				for(int k = 0; k < 3; k++)
					owner.xlats[j][k] += diffxlat[k];
				Transformation temp = new Transformation(owner.axes[j], owner.angles[j], owner.xlats[j]);
				temp.rotMatrix = diffRotMat.applyTo(temp.rotMatrix);
				owner.axes[j] = temp.rotMatrix.getAxis().toArray();
				owner.angles[j] = temp.rotMatrix.getAngle();
			}
		}
		
		return logProposalRatio;
	}

}
