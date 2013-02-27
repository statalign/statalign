package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.linear.ArrayRealVector;

import statalign.model.ext.plugins.StructAlign;
import statalign.model.ext.plugins.structalign.Transformation;

public class LibraryMove extends RotationOrTranslationMove {

	public LibraryMove (StructAlign s, String n) {
		owner = s; // As an McmcModule object
		structAlign = s;
		// Since Java can't handle templates, we use this
		// variable to avoid multiple casts throughout the code.
		
		name = n;
	}

	public double proposal(Object externalState) {
		
		structAlign.rotProp.sd = proposalWidthControlVariable;
		Transformation oldSub = new Transformation(structAlign.axes[index], structAlign.angles[index], structAlign.xlats[index]);
		// transformation should be relative to reference protein
		oldSub.xlat = oldSub.xlat.subtract(new ArrayRealVector(structAlign.xlats[0]));
		Transformation libProp = structAlign.rotProp.propose(index);
		structAlign.axes[index] = libProp.rotMatrix.getAxis().toArray();
		structAlign.angles[index] = libProp.rotMatrix.getAngle();
		structAlign.xlats[index] = libProp.xlat.toArray();

		// library density 
		double logProposalRatio = structAlign.rotProp.libraryLogDensity(index, oldSub) - 
				structAlign.rotProp.libraryLogDensity(index, libProp);

		// proposed translation is relative to reference protein
		for(int i = 0; i < 3; i++)
			structAlign.xlats[index][i] += structAlign.xlats[0][i];			

		// calculate 'difference' between proposed and current transformations
		double[] diffxlat = new double[3];
		for(int i = 0; i < 3; i++)
			diffxlat[i] = structAlign.xlats[index][i] - oldxlats[index][i];
		Rotation diffRotMat = libProp.rotMatrix.applyTo(oldSub.rotMatrix.revert());

		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			if(j != index){
				for(int k = 0; k < 3; k++)
					structAlign.xlats[j][k] += diffxlat[k];
				Transformation temp = new Transformation(structAlign.axes[j], structAlign.angles[j], structAlign.xlats[j]);
				temp.rotMatrix = diffRotMat.applyTo(temp.rotMatrix);
				structAlign.axes[j] = temp.rotMatrix.getAxis().toArray();
				structAlign.angles[j] = temp.rotMatrix.getAngle();
			}
		}
		
		return logProposalRatio;
	}

}
