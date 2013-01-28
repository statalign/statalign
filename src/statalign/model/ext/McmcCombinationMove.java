package statalign.model.ext;

import java.util.List;

public abstract class McmcCombinationMove extends McmcMove {

	protected List<McmcMove> mcmcMoves;
	protected List<Integer> mcmcMoveWeights;
	
	protected double oldll;
	
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
	public abstract void updateLikelihood(Object externalState);
	
	public void restoreState(Object externalState) {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.restoreState(externalState);
		}
	}
}
