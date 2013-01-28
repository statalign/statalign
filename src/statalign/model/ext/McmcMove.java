package statalign.model.ext;

import statalign.model.ext.ModelExtension;

public abstract class McmcMove {

	ModelExtension owner;
	public String name;
	//public T value;
	public int proposalCount;
	public int acceptanceCount;
	public boolean lastMoveAccepted;
	public double proposalWidthControlVariable;
	public boolean autoTune;
	
	public double acceptanceRate() {
		return (double) acceptanceCount / (double) proposalCount;
	}
	
	public abstract void copyState();
	public abstract double proposal(Object externalState); // Returns logProposalRatio
	public abstract void updateLikelihood(Object externalState); 
	public abstract void restoreState();
	
	public void move(Object externalState) {
		
		proposalCount++;
		copyState();
		double logProposalRatio = proposal(externalState);
		updateLikelihood(externalState);
		if(owner.isParamChangeAccepted(logProposalRatio)) {
			acceptanceCount++;
			lastMoveAccepted = true;
		}
		else {
			lastMoveAccepted = false;
			restoreState();
		}
	}
	
	 
	
}
