package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.linear.ArrayRealVector;

import statalign.model.ext.plugins.StructAlign;
import statalign.model.ext.plugins.structalign.Transformation;

public class LibraryMove extends RotationOrTranslationMove {

	StructAlign owner;
	
	public LibraryMove (StructAlign s) {
		owner = s;
		omit = 1;
	}

	public double proposal(Object externalState) {
		Transformation old = new Transformation(owner.axes[ind], owner.angles[ind], owner.xlats[ind]);
		// transformation should be relative to reference protein
		old.xlat = old.xlat.subtract(new ArrayRealVector(owner.xlats[0]));
		Transformation prop = owner.rotProp.propose(ind);
		owner.axes[ind] = prop.rotMatrix.getAxis().toArray();
		owner.angles[ind] = prop.rotMatrix.getAngle();
		owner.xlats[ind] = prop.xlat.toArray();

		// proposed translation is relative to reference protein
		for(int i = 0; i < 3; i++)
			owner.xlats[ind][i] += owner.xlats[0][i];
		
		return (owner.rotProp.libraryLogDensity(ind, old) - 
				  owner.rotProp.libraryLogDensity(ind, prop));
	}

}
