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
	double P; // for hmm3
	double proposalParamMultiplier = 1.0;
	double[] realParams;
	
	public boolean autoTunable = true;
	private boolean useModelExtInProposal = false;
	private boolean useModextEm = false, useModextUpp = false;
	
	double oldll; 
	
	final static double LEAFCOUNT_POWER = 1.0; // Original
	//final static double LEAFCOUNT_POWER = -2.0;
	final static double SELTRLEVPROB[] = { 0.9, 0.6, 0.4, 0.2, 0 };
	
	public double minAcceptance = 0.05; // keep tuning till we get to this
	public static final double MIN_WINDOW_MULTIPLIER = 0.5;
	public static final double MAX_WINDOW_MULTIPLIER = 1.5;
	
	public AlignmentMove (McmcModule m, double _P, String n) {
		owner = m;
		name = n;		
		P = _P;
		spanMultiplier = 0.9;
		autoTune = false;
	}
	public void setProposalParamMultiplier(double p) {
		proposalParamMultiplier = p;
	}
	public void useModelExtInProposal() {
		useModelExtInProposal = true;
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
		if (!useModelExtInProposal) {
			useModextEm = Utils.USE_MODEXT_EM; 
			useModextUpp = Utils.USE_MODEXT_UPP; 
			Utils.USE_MODEXT_EM = false;			
			Utils.USE_MODEXT_UPP = false;
		}		
		if (proposalWidthControlVariable < MIN_WINDOW_MULTIPLIER) {		
			proposalWidthControlVariable = MIN_WINDOW_MULTIPLIER; 
		}		
		if (proposalWidthControlVariable > MAX_WINDOW_MULTIPLIER) {		
			proposalWidthControlVariable = MAX_WINDOW_MULTIPLIER;
		}
		
		Utils.WINDOW_MULTIPLIER = proposalWidthControlVariable;
		
		oldll = owner.curLogLike;
		
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
				
		tree.hmm3.updateParam(new double[]{P});				
		selectedRoot.recursivelyUpdateHmm3Matrices();
		
		if (proposalParamMultiplier != 1.0) {
			realParams = tree.hmm2.params.clone();
			if (selectedRoot != tree.root) {
				selectedRoot.updateHmm2Matrix(new double[] {realParams[0], 
						 proposalParamMultiplier*realParams[1], 
						 proposalParamMultiplier*realParams[2]});
			}
		}
	}
	public double proposal(Object externalState) {
		
		((CoreMcmcModule) owner).getModelExtMan().beforeAlignChange(tree, selectedRoot);
		double logProposalRatio = selectedRoot.selectAndResampleAlignment();
		if (proposalParamMultiplier != 1.0 && selectedRoot != tree.root) {
			selectedRoot.updateHmm2Matrix(realParams);
			selectedRoot.calcOrphan();
			selectedRoot.parent.calcFelsen();
			selectedRoot.parent.calcOrphan();
			selectedRoot.parent.calcIndelLogLike();
			selectedRoot.calcAllUp();
		}
		// NB need to reset them here, because they're used in the likelihood computation as well
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
		selectedRoot.calcOrphan();
		if (selectedRoot != tree.root) {
			selectedRoot.parent.calcFelsen();
			selectedRoot.parent.calcOrphan();
			selectedRoot.parent.calcIndelLogLike();
			selectedRoot.calcAllUp();
		}
//		 tree.root.calcFelsenRecursively();
//         tree.root.calcOrphanRecursively();
//         tree.root.calcIndelLogLikeRecursively();
		
//		selectedRoot.calcFelsenRecursively();
//		selectedRoot.calcOrphanRecursively();
//		selectedRoot.calcIndelLogLikeRecursively();
         if (Utils.USE_UPPER) {
         	//owner.root.calcFelsenRecursively();
         	//tree.root.calcUpperRecursively();
        	selectedRoot.calcUpperFromRoot();
         }   
	}
	
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterAlignChange(tree, selectedRoot,lastMoveAccepted);
		if (lastMoveAccepted && (owner.curLogLike == oldll)) acceptanceCount--;
				
		if (!useModelExtInProposal) {			
			Utils.USE_MODEXT_EM = useModextEm;			
			Utils.USE_MODEXT_UPP = useModextUpp;
			if (lastMoveAccepted) {
				selectedRoot.updateAlignedRecursively();
				selectedRoot.updateAlignedParent();
			}
		}	
	}
	
	@Override
	public void afterFirstHalfBurnin() {
		if (autoTunable) {
			autoTune = true;
		}
	}
	 
	
}
