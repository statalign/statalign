package statalign.model.ext;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Tree;

public class McmcCombinationMove extends McmcMove {

	protected List<McmcMove> mcmcMoves;
	protected List<Integer> mcmcMoveWeights;
	
	Tree tree;
	
	public McmcCombinationMove (ArrayList<McmcMove> mcmcMoves) {
		if (mcmcMoves.size() < 2) {
			throw new IllegalArgumentException("McmcCombinationMove must contain at least two McmcMove objects");
		}
		
		name = mcmcMoves.get(0).name;
		for (int i=1; i<mcmcMoves.size(); i++) {
			name += "_"+mcmcMoves.get(i).name;
		}
	}
	public void copyState(Object externalState) {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.copyState(externalState);
		}
	}
	public double proposal(Object externalState) {
		double logProposalRatio = 0;
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.copyState(externalState);
			logProposalRatio -= mcmcMove.logPriorDensity(externalState);
			logProposalRatio += mcmcMove.proposal(externalState); 
			logProposalRatio += mcmcMove.logPriorDensity(externalState);
		}
		return logProposalRatio;
	}
	
	public double logPriorDensity(Object externalState) {
		double logPriorDensity = 0;
		for (McmcMove mcmcMove : mcmcMoves) {
			logPriorDensity += mcmcMove.logPriorDensity(externalState);
		}
		return logPriorDensity;
	}
	public void updateLikelihood(Object externalState) {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.updateLikelihood(externalState);
		}
	}
	
	public void restoreState(Object externalState) {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.restoreState(externalState);
		}
	}
}
