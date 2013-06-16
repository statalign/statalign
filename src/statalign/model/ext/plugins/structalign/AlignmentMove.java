package statalign.model.ext.plugins.structalign;

import java.util.ArrayList;

import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.Tree;
import statalign.mcmc.McmcMove;
import statalign.model.ext.plugins.StructAlign;

public class AlignmentMove extends McmcMove {

	Tree tree;
	public StructAlignMoveParams moveParams = new StructAlignMoveParams();
	
	double[][] oldaxes = null;
	double[] oldangles = null;
	double[][] oldxlats = null;
	double[][][] oldrots = null;
	double oldll;
	
	double heat = 1.0;
	double secondHeat = 1.0;
	
	Vertex subtreeRoot;
	int nLeaves;
	ArrayList<Integer> subtreeLeaves;
	int index;
	
	public AlignmentMove (StructAlign s, String n) {
		owner = s;
		name = n;
		autoTune = false; 
		// This move gets autoTune'd via the core AlignmentMove
	}
	
	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		else {
			throw new IllegalArgumentException("AlignmentMove.copyState must take an argument of type Tree.");
		}
		subtreeRoot = Funcs.sampleVertex(tree);
		nLeaves = ((StructAlign) owner).coords.length;
		subtreeLeaves = Subtree.getSubtreeLeaves(tree, subtreeRoot, nLeaves);
		index = subtreeLeaves.get(Utils.generator.nextInt(subtreeLeaves.size()));
		((StructAlign) owner).oldAlign = ((StructAlign) owner).curAlign;
		oldll = owner.getLogLike();
	}

	public double proposal(Object externalState) {
		double logProposalRatio = subtreeRoot.realignToParent(heat);
		((StructAlign) owner).curAlign = tree.getState().getLeafAlign();
		return logProposalRatio;
	}
	
	public double logPriorDensity(Object externalState) {
		return 0; // Uniform prior
	}

	public void updateLikelihood(Object externalState) {		
		owner.setLogLike( ((StructAlign) owner).calcAllColumnContrib() );
	}
	public void restoreState(Object externalState) {
		subtreeRoot.alignRestore(); 
		
		((StructAlign) owner).curAlign = ((StructAlign) owner).oldAlign;
		owner.setLogLike(oldll);
	}
	public void afterMove(Object externalState) {
		if (Utils.DEBUG) {
			tree.root.calcFelsenRecursively(); 
			tree.root.calcOrphanRecursively(); 
			tree.root.calcIndelLogLikeRecursively(); 
	        if (Utils.USE_UPPER) {
	        	//owner.root.calcFelsenRecursively();
	        	tree.root.calcUpperRecursively();
	        }   
	        //tree.root.recomputeCheckLogLike();
		}
		//if (lastMoveAccepted && (owner.getLogLike() == oldll)) acceptanceCount--;
	}
	@Override
	public void afterFirstHalfBurnin() {
		heat = secondHeat;
	}
}
