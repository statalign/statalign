package statalign.mcmc;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Utils;
import static statalign.utils.Reversed.reversed;

public class McmcCombinationMove extends McmcMove {

	protected List<McmcMove> mcmcMoves = new ArrayList<McmcMove>();
	
	public McmcCombinationMove (ArrayList<McmcMove> moves) {
		if (moves.size() < 2) {
			throw new IllegalArgumentException("McmcCombinationMove must contain at least two McmcMove objects");
		}		
		name = moves.get(0).name;
		mcmcMoves.add(moves.get(0));
		for (int i=1; i<moves.size(); i++) {
			name += "_"+moves.get(i).name;
			mcmcMoves.add(moves.get(i));
		}
	}
	public void copyState(Object externalState) {
		if (Utils.DEBUG) System.out.println(name+" "+proposalWidthControlVariable);
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.inCombination = true;
			mcmcMove.copyState(externalState);
			mcmcMove.inCombination = false;
		}
	}
	public double proposal(Object externalState) {
		double logProposalRatio = 0;
		for (McmcMove mcmcMove : mcmcMoves) {
			if (Utils.DEBUG) System.out.print(mcmcMove.name+" "+mcmcMove.proposalWidthControlVariable+" ");
			mcmcMove.inCombination = true;			
			double oldProposalWidthControlVariable = mcmcMove.proposalWidthControlVariable;
			mcmcMove.proposalWidthControlVariable *= proposalWidthControlVariable;
			//mcmcMove.proposalCount++;
			logProposalRatio += mcmcMove.proposal(externalState); 
			mcmcMove.proposalWidthControlVariable = oldProposalWidthControlVariable;
			mcmcMove.inCombination = false;
			// In combination moves, the counts of the sub-moves are not incremented
		}
		if (Utils.DEBUG) System.out.println();
		// We implicitly assume here that the order of the proposals
		// does not matter.
		return logProposalRatio;
	}
	
	@Override
	public McmcModule getOwner() { 
		return mcmcMoves.get(0).getOwner();
		// NB it doesn't matter which McmcMove we use here, since it is
		// only used to call back to the Mcmc object running
		// all of them, but since the McmcCombinationMove does not necessarily
		// have a unique ModelExtension as its owner, we must use one of
		// its McmcMove objects to do the callback.
	}
			
	public double logPriorDensity(Object externalState) {
		double logPriorDensity = 0;
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.inCombination = true;
			logPriorDensity += mcmcMove.logPriorDensity(externalState);
			mcmcMove.inCombination = false;
		}
		return logPriorDensity;
	}
	public void updateLikelihood(Object externalState) {
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.inCombination = true;
			mcmcMove.updateLikelihood(externalState);
			mcmcMove.inCombination = false;
		}
	}
	public void afterMove(Object externalState) {
		if (Utils.DEBUG) System.out.println((this.lastMoveAccepted ? "Accepted (" : "Rejected (") + this.acceptanceRate() + ")");		
		for (McmcMove mcmcMove : mcmcMoves) {
			mcmcMove.inCombination = true;			
			mcmcMove.afterMove(externalState);			
			mcmcMove.inCombination = false;
		}
	}
	public void restoreState(Object externalState) {
		for (McmcMove mcmcMove : reversed(mcmcMoves)) {
			mcmcMove.inCombination = true;
			mcmcMove.restoreState(externalState);
			mcmcMove.inCombination = false;
		}
	}
}
