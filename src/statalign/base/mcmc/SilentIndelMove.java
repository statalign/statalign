package statalign.base.mcmc;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;

public class SilentIndelMove extends McmcMove {

	Tree tree = null;
	Vertex v;	
	double[] weights;
	final static double LEAFCOUNT_POWER = 1.0; // Original
	final static double SELTRLEVPROB[] = { 0.9, 0.6, 0.4, 0.2, 0 };
	
	boolean didInsertion = false;
	boolean validProposal = false;
	
	public SilentIndelMove (McmcModule m, String n) {
		owner = m;
		name = n;				
	}

	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else throw new IllegalArgumentException("SilentIndelMove.copyState must take an argument of type Tree.");
		weights = new double[tree.vertex.length];		
		tree.countLeaves();
		for (int i = 0; i < weights.length; i++) {
			weights[i] = Math.pow(tree.vertex[i].leafCount, LEAFCOUNT_POWER);
		}
		int k = 0;
		v = null;
		System.out.println("Root = "+tree.root);
		while (v==null || v.leafCount==1) { // misbehaves with saveFiveWay -- may be no uncle
		//while (v==null || v == tree.root || v.parent == tree.root || v.leafCount==1) {
			k = Utils.weightedChoose(weights, null);			
			v = tree.vertex[k];
			System.out.println(k+" "+v.parent);
			v.selectSubtree(SELTRLEVPROB, 0);
		}
		if (Utils.DEBUG) System.out.println("Vertex = "+k);
		v.saveData();
		((CoreMcmcModule) owner).getModelExtMan().beforeAlignChange(tree, v);				
	}
	public double proposal(Object externalState) {
				
		
		double logProposalRatio = 0;
		if (Utils.generator.nextDouble() < Utils.SILENT_INSERT_PROB) {
			if (Utils.DEBUG) System.out.println("Inserting silent indel.");
			didInsertion = true;
			logProposalRatio = v.insertSilentIndel();
		}
		else {
			if (Utils.DEBUG) System.out.println("Excising silent indel.");
			didInsertion = false;
			logProposalRatio = v.exciseSilentIndel(); 
		}
		if (logProposalRatio == Double.NEGATIVE_INFINITY) {
			validProposal = false;
		}
		else validProposal = true;
			
		return logProposalRatio;
	}
	public double logPriorDensity(Object externalState) {
		return 0.0;
	}
	public void updateLikelihood(Object externalState) {
		if (validProposal) {
			v.calcFelsen();
            v.calcOrphan();
            v.calcIndelLogLike();
			v.calcAllUp();
			owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeAlignChange(tree, v));
		}
	}
	public void restoreState(Object externalState) {
//		v.restoreFiveWay(false); // currently this is the way we're saving the old state
		if (validProposal) {
	 		if (didInsertion) v.undoInsertSilentIndel();
			else v.undoExciseSilentIndel();	 		
		}
		//v.restoreData();
		v.calcFelsen();
        v.calcOrphan();
        v.calcIndelLogLike();
		v.calcAllUp();
		// this could be problematic if the proposal failed
	}
	
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterAlignChange(tree, v,lastMoveAccepted);
	}
	
	 
	
}
