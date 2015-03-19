package statalign.model.ext.plugins.structalign;

import java.util.ArrayList;

import org.apache.commons.math3.util.MathArrays;

import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.Tree;
import statalign.mcmc.McmcMove;
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
	AlignmentMove alignmentMove;
	int nLeaves;
	ArrayList<Integer> subtreeLeaves;
	int nStrucInSubtree;
	int index;

	public void shareSubtreeRoot(AlignmentMove m) {
		alignmentMove = m;
		System.out.println("Move '"+name+"' now sharing subtree root with move '"+m.name+"'");
	}
	public Vertex getSharedSubtreeRoot() {
		return alignmentMove.getSubtreeRoot();
	}
	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		else {
			throw new IllegalArgumentException("RotationOrTranslationMove.copyState must take an argument of type Tree.");
		}
		if (Utils.DEBUG) {
			System.out.print(inCombination?"In combination\n":"");
			System.out.println("Shared root = "+getSharedSubtreeRoot());
		}
		if (inCombination) subtreeRoot = getSharedSubtreeRoot();
		else subtreeRoot = Funcs.sampleVertex(tree);		
		nLeaves = structAlign.coords.length;
		subtreeLeaves = Subtree.getSubtreeLeaves(tree, subtreeRoot, nLeaves);
		index = subtreeLeaves.get(Utils.generator.nextInt(subtreeLeaves.size()));
		int refInd = 0;
		//while (structAlign.axes[refInd]==null) ++refInd;
		while (structAlign.coords[refInd]==null) ++refInd;
		oldaxes = new double[structAlign.axes.length][structAlign.axes[refInd].length];
		oldangles = new double[structAlign.angles.length];
		oldxlats = new double[structAlign.xlats.length][structAlign.xlats.length];
		oldrots = new double[structAlign.rotCoords.length][structAlign.rotCoords[refInd].length][structAlign.rotCoords[refInd][refInd].length];
		nStrucInSubtree = 0;
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			//if (structAlign.axes[j] == null) continue;
			if (structAlign.coords[j] == null) continue;
			nStrucInSubtree++;
			oldaxes[j] = MathArrays.copyOf(structAlign.axes[j]);
			oldangles[j] = structAlign.angles[j];
			oldxlats[j] = MathArrays.copyOf(structAlign.xlats[j]);
			oldrots[j] = structAlign.rotCoords[j];
		}
		
		oldll = owner.getLogLike();
	}

	public double proposal(Object externalState) {
//		if (nStrucInSubtree==0) return Double.NEGATIVE_INFINITY;
//		else return 0;
		//System.out.println(name+" "+proposalWidthControlVariable);
		return 0;
	}
	
	public double logPriorDensity(Object externalState) {
		return 0; // Uniform prior
	}

	public void updateLikelihood(Object externalState) {
		for(int i = 0; i < subtreeLeaves.size(); i++){
			int j = subtreeLeaves.get(i);
			if (structAlign.coords[j] == null) continue;
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
	public void afterMove(Object externalState) {			
	}
}
