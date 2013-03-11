package statalign.base.mcmc;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;

public class AlignmentMove extends McmcMove {

	Tree tree = null;
	Vertex selectedRoot;
	double[] weights;
	final static double LEAFCOUNT_POWER = 1.0;
	final static double SELTRLEVPROB[] = { 0.9, 0.6, 0.4, 0.2, 0 };
	
	public static final double MIN_WINDOW_MULTIPLIER = 0.1;
	public static final double MAX_WINDOW_MULTIPLIER = 1.0;
	
	public AlignmentMove (McmcModule m, String n) {
		owner = m;
		name = n;		
		spanMultiplier = 0.9;
	}

	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("AlignmentMove.copyState must take an argument of type Tree.");
		}
		if (proposalWidthControlVariable > MIN_WINDOW_MULTIPLIER &&
				proposalWidthControlVariable < MAX_WINDOW_MULTIPLIER) {
			Utils.WINDOW_MULTIPLIER = proposalWidthControlVariable;
		}
	}
	public double proposal(Object externalState) {
		weights = new double[tree.vertex.length];
		for (int i = 0; i < tree.vertex.length; i++) {
			tree.vertex[i].selected = false;
		}
		tree.countLeaves();
		for (int i = 0; i < weights.length; i++) {
			weights[i] = Math.pow(tree.vertex[i].leafCount, LEAFCOUNT_POWER);
		}
		int k = Utils.weightedChoose(weights, null);
		selectedRoot = tree.vertex[k];
		selectedRoot.selectSubtree(SELTRLEVPROB, 0);
		((CoreMcmcModule) owner).getModelExtMan().beforeAlignChange(tree, selectedRoot);
		double logProposalRatio = selectedRoot.selectAndResampleAlignment();
		return logProposalRatio;
	}
	public double logPriorDensity(Object externalState) {
		return 0.0;
	}
	public void updateLikelihood(Object externalState) {
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeAlignChange(tree, selectedRoot));
	}
	public void restoreState(Object externalState) {
		selectedRoot.alignRestore();
	}
	
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterAlignChange(tree, selectedRoot,lastMoveAccepted);
	}
	
	 
	
}
