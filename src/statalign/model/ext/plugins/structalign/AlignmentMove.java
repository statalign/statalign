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
	
	boolean autoTunable = true;	
	
	public Vertex subtreeRoot;
	int nLeaves;
	ArrayList<Integer> subtreeLeaves;
	int index;
	
	public double minAcceptance = 0.05; // keep tuning till we get to this
	public static final double MIN_WINDOW_MULTIPLIER = 0.5;
	public static final double MAX_WINDOW_MULTIPLIER = 1.5;
	
	public AlignmentMove (StructAlign s, String n) {
		owner = s;
		name = n;
		autoTune = false; 
		proposalWidthControlVariable = 1.0;		
		// This move gets autoTune'd via the core AlignmentMove
	}

	public Vertex getSubtreeRoot() { 
		return subtreeRoot;
	}
	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		else {
			throw new IllegalArgumentException("AlignmentMove.copyState must take an argument of type Tree.");
		}
		if (proposalWidthControlVariable < MIN_WINDOW_MULTIPLIER) {		
			proposalWidthControlVariable = MIN_WINDOW_MULTIPLIER; 
		}		
		if (proposalWidthControlVariable > MAX_WINDOW_MULTIPLIER) {		
			proposalWidthControlVariable = MAX_WINDOW_MULTIPLIER;
		}
		
		Utils.WINDOW_MULTIPLIER = proposalWidthControlVariable;
		
		subtreeRoot = Funcs.sampleVertex(tree);		
		if (Utils.DEBUG) System.out.println("subtreeRoot = "+subtreeRoot.index);
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
		if (Utils.DEBUG && Utils.USE_MODEXT_EM) tree.root.updateAlignedRecursivelyWithCheck();
		if (Utils.DEBUG) tree.root.recomputeCheckLogLike();
		subtreeRoot.alignRestore();
//		subtreeRoot.calcOrphan();
//		subtreeRoot.calcAllUp();
		//System.out.println("oldll = "+oldll+", newll = "+owner.getLogLike());
		if (Utils.DEBUG) tree.root.recomputeCheckLogLike();
		if (Utils.USE_MODEXT_EM) subtreeRoot.updateAlignedParentInWindow();
		if (Utils.DEBUG && Utils.USE_MODEXT_EM) tree.root.updateAlignedRecursivelyInWindowWithCheck();
		if (Utils.DEBUG && Utils.USE_MODEXT_EM) tree.root.updateAlignedRecursivelyWithCheck();
		
		((StructAlign) owner).curAlign = ((StructAlign) owner).oldAlign;
		owner.setLogLike(oldll);
	}
	public void afterMove(Object externalState) {
		if (lastMoveAccepted && lastLogProposalRatio == 0) acceptanceCount--;
		
		if (Utils.DEBUG) {
//			tree.root.calcFelsenRecursively(); 
//			tree.root.calcOrphanRecursively(); 
//			tree.root.calcIndelLogLikeRecursively(); 
//	        if (Utils.USE_UPPER) {
//	        	//owner.root.calcFelsenRecursively();
//	        	subtreeRoot.calcUpperFromRoot();
//	        	//tree.root.calcUpperRecursively();
//	        }   
	        //tree.root.recomputeCheckLogLike();
		}
//		if (lastMoveAccepted && subtreeRoot != tree.root) {						
//			subtreeRoot.updateAlignedParent();
//		}
		//if (lastMoveAccepted && (owner.getLogLike() == oldll)) acceptanceCount--;
	}
	@Override
	public void afterFirstHalfBurnin() {
		heat = secondHeat;
		//autoTune = autoTunable;
	}
}
