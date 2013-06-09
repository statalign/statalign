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
	//final static double LEAFCOUNT_POWER = 2.0; 
	final static double LEAFCOUNT_POWER = 1.0; // Original
	
	boolean didInsertion = false;
	boolean validProposal = false;
	
	public SilentIndelMove (McmcModule m, String n) {
		owner = m;
		name = n;				
		autoTune = false;
	}

	public void copyState(Object externalState) {
		Utils.DEBUG = true;
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else throw new IllegalArgumentException("SilentIndelMove.copyState must take an argument of type Tree.");
		weights = new double[tree.vertex.length];		
		tree.countLeaves();
		tree.countSilentIndels();
		for (int i = 0; i < weights.length; i++) {
			//weights[i] = Math.pow(tree.vertex[i].leafCount, LEAFCOUNT_POWER);
			//weights[i] = (tree.vertex[i].nSilentIndels > 0) ? 1 : 0;
			weights[i] = (tree.vertex[i].nSilentIndels > 1) ? 1 : 0;
		}
		int k = 0;
		v = null;
		if (Utils.DEBUG) System.out.println("Root = "+tree.root);
		while (v==null || v.leafCount==1) { // misbehaves with saveFiveWay -- may be no uncle
		//while (v==null || v == tree.root || v.parent == tree.root || v.leafCount==1) {
			k = Utils.weightedChoose(weights, null);			
			v = tree.vertex[k];
			if (Utils.DEBUG) System.out.println(k+" "+v.parent);			
		}
		if (Utils.DEBUG) System.out.println("Vertex = "+k);
		v.saveData();
		((CoreMcmcModule) owner).getModelExtMan().beforeAlignChange(tree, v);				
	}
	public double proposal(Object externalState) {
				
		if (Utils.DEBUG) tree.root.printToScreenAlignment(0,0,true);
		
		double logProposalRatio = 0;		
		if (Utils.generator.nextDouble() < Utils.SILENT_INSERT_PROB) {
			if (Utils.DEBUG) System.out.println("Inserting silent indel.");
			didInsertion = true;			
			//logProposalRatio = v.insertSilentIndel();
			logProposalRatio = v.modifySilentIndel(didInsertion);
		}
		else {
			if (Utils.DEBUG) System.out.println("Excising silent indel.");
			didInsertion = false;
			//logProposalRatio = v.exciseSilentIndel(); 
			logProposalRatio = v.modifySilentIndel(didInsertion);
		}
		if (logProposalRatio == Double.NEGATIVE_INFINITY) {
			validProposal = false;
		}
		else validProposal = true;

		if (Utils.DEBUG) tree.root.printToScreenAlignment(0,0,true);
		
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
		Utils.DEBUG = false;
		((CoreMcmcModule) owner).getModelExtMan().afterAlignChange(tree, v,lastMoveAccepted);
	}
	
	 
	
}
