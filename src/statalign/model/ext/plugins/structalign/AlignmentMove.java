package statalign.model.ext.plugins.structalign;

import java.util.ArrayList;

import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.Tree;
import statalign.model.ext.plugins.StructAlign;

public class AlignmentMove extends StructAlignMcmcMove {

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
	
	public AlignmentMove (StructAlign s, String n) {
		owner = s;
		name = n;
	}
	
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
		owner.oldAlign = owner.curAlign;
		oldll = owner.curLogLike;
	}

	public double proposal(Object externalState) {
		double logProposalRatio = subtreeRoot.realignToParent();
		owner.curAlign = tree.getState().getLeafAlign();
		return logProposalRatio;
	}
	
	public double logPriorDensity(Object externalState) {
		return 0; // Uniform prior
	}

	public void updateLikelihood(Object externalState) {		
		owner.curLogLike = owner.calcAllColumnContrib();
	}
	public void restoreState(Object externalState) {
		subtreeRoot.alignRestore();
		owner.curAlign = owner.oldAlign;
		owner.curLogLike = oldll;
	}
}
