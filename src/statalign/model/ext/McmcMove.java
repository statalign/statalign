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
	public abstract double proposal(); // Returns logProposalRatio
	public abstract void updateLikelihood(); 
	public abstract void restoreState();
	
	public void move() {
		
		proposalCount++;
		copyState();
		double logProposalRatio = proposal();
		updateLikelihood();
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
