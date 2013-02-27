package statalign.model.ext.plugins.structalign;

import java.util.ArrayList;

import org.apache.commons.math3.util.MathArrays;

import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.Tree;
import statalign.base.mcmc.McmcMove;
import statalign.model.ext.plugins.StructAlign;

public abstract class RotationOrTranslationMove extends McmcMove {

	Tree tree;
	protected StructAlign structAlign;
	public StructAlignMoveParams moveParams = new StructAlignMoveParams(); 

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
		nLeaves = structAlign.coords.length;
		subtreeLeaves = Subtree.getSubtreeLeaves(tree, subtreeRoot, nLeaves);
		index = subtreeLeaves.get(Utils.generator.nextInt(subtreeLeaves.size()));
		oldaxes = new double[structAlign.axes.length][structAlign.axes[0].length];
		oldangles = new double[structAlign.angles.length];
		oldxlats = new double[structAlign.xlats.length][structAlign.xlats.length];
		oldrots = new double[structAlign.rotCoords.length][structAlign.rotCoords[0].length][structAlign.rotCoords[0][0].length];
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			oldaxes[j] = MathArrays.copyOf(structAlign.axes[j]);
			oldangles[j] = structAlign.angles[j];
			oldxlats[j] = MathArrays.copyOf(structAlign.xlats[j]);
			oldrots[j] = structAlign.rotCoords[j];
		}
		
		oldll = owner.getLogLike();
	}

	public abstract double proposal(Object externalState);
	
	public double logPriorDensity(Object externalState) {
		return 0; // Uniform prior
	}

	public void updateLikelihood(Object externalState) {
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			structAlign.rotCoords[j] = null;	// so that calcRotation creates new array
			structAlign.calcRotation(j);
		}
		owner.setLogLike( structAlign.calcAllColumnContrib() );
	}
	public void restoreState(Object externalState) {
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			structAlign.axes[j] = oldaxes[j];
			structAlign.angles[j] = oldangles[j];
			structAlign.xlats[j] = oldxlats[j];
			structAlign.rotCoords[j] = oldrots[j];
		}
		owner.setLogLike(oldll);
	}
}
