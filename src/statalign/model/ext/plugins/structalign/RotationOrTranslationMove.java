package statalign.model.ext.plugins.structalign;

import java.util.ArrayList;

import org.apache.commons.math3.util.MathArrays;

import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.Tree;

public abstract class RotationOrTranslationMove extends StructAlignMcmcMove {

	Tree tree;
	
	double[][] oldaxes = null;
	double[] oldangles = null;
	double[][] oldxlats = null;
	double[][][] oldrots = null;
	double oldll;
	
	Vertex subtreeRoot;
	int nLeaves;
	ArrayList<Integer> subtreeLeaves;
	int index;

	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		else {
			throw new IllegalArgumentException("RotationOrTranslationMove.copyState must take an argument of type Tree.");
		}
		subtreeRoot = Funcs.sampleVertex(tree);
		nLeaves = owner.coords.length;
		subtreeLeaves = Subtree.getSubtreeLeaves(tree, subtreeRoot, nLeaves);
		index = subtreeLeaves.get(Utils.generator.nextInt(subtreeLeaves.size()));
		oldaxes = new double[owner.axes.length][owner.axes[0].length];
		oldangles = new double[owner.angles.length];
		oldxlats = new double[owner.xlats.length][owner.xlats.length];
		oldrots = new double[owner.rotCoords.length][owner.rotCoords[0].length][owner.rotCoords[0][0].length];
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			oldaxes[j] = MathArrays.copyOf(owner.axes[j]);
			oldangles[j] = owner.angles[j];
			oldxlats[j] = MathArrays.copyOf(owner.xlats[j]);
			oldrots[j] = owner.rotCoords[j];
		}
		
		oldll = owner.curLogLike;
	}

	public abstract double proposal(Object externalState);
	
	public double logPriorDensity(Object externalState) {
		return 0; // Uniform prior
	}

	public void updateLikelihood(Object externalState) {
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			owner.rotCoords[j] = null;	// so that calcRotation creates new array
			owner.calcRotation(j);
		}
		owner.curLogLike = owner.calcAllColumnContrib();
	}
	public void restoreState(Object externalState) {
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			owner.axes[j] = oldaxes[j];
			owner.angles[j] = oldangles[j];
			owner.xlats[j] = oldxlats[j];
			owner.rotCoords[j] = oldrots[j];
		}
		owner.curLogLike = oldll;
	}
}
