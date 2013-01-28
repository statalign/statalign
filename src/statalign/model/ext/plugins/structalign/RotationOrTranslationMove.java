package statalign.model.ext.plugins.structalign;

import org.apache.commons.math3.util.MathArrays;

import statalign.base.Utils;
import statalign.model.ext.McmcMove;
import statalign.model.ext.plugins.StructAlign;

public abstract class RotationOrTranslationMove extends McmcMove {

	StructAlign owner;
	
	int omit; // number of structures to be omitted
	int ind; // index of structure to be rotated or translated
	
	double[] oldax;
	double oldang;
	double[] oldxlat;
	double[][] oldrots;
	double oldll;
	
	public void copyState() {
		ind = Utils.generator.nextInt(owner.coords.length - omit) + omit;
		
		oldax = MathArrays.copyOf(owner.axes[ind]);
		oldang = owner.angles[ind];
		oldxlat = MathArrays.copyOf(owner.xlats[ind]);
		oldrots = owner.rotCoords[ind];
		owner.rotCoords[ind] = null;	// so that calcRotation creates new array
		oldll = owner.curLogLike;
	}

	public abstract double proposal(Object externalState);
	
	public void updateLikelihood(Object externalState) {
		owner.calcRotation(ind);
		owner.curLogLike = owner.calcAllColumnContrib();
	}
	public void restoreState() {
		owner.axes[ind] = oldax;
		owner.angles[ind] = oldang;
		owner.xlats[ind] = oldxlat;
		owner.rotCoords[ind] = oldrots;
		owner.curLogLike = oldll;
	}
}
