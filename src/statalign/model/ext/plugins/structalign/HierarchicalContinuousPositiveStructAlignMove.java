package statalign.model.ext.plugins.structalign;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Utils;
import statalign.mcmc.GammaPrior;
import statalign.mcmc.LogNormalPrior;
import statalign.mcmc.McmcMove;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class HierarchicalContinuousPositiveStructAlignMove extends ContinuousPositiveStructAlignMove {

	private List<ContinuousPositiveStructAlignMove> children = new ArrayList<ContinuousPositiveStructAlignMove>();
	public PriorDistribution<Double> hierarchicalPrior;
	private boolean onlySampleIfAtLeastTwoChildrenNotFixed = false;
	
	boolean allowSpikeSelection = false;
	public void disallowSpikeSelection() { allowSpikeSelection = false; }
	public void allowSpikeSelection() { allowSpikeSelection = true; }
	double[] spikeParams;
	
	public HierarchicalContinuousPositiveStructAlignMove (StructAlign s, 
			ParameterInterface p,  PriorDistribution<Double> pr,
			ProposalDistribution<Double> prop, String n) {
		super(s,p,pr,prop,n);
		//hierarchicalPrior = new GammaPrior(((StructAlign) owner).nu,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		hierarchicalPrior = new LogNormalPrior(Math.log(((StructAlign) owner).sigma2Hier),((StructAlign) owner).nu);
		//hierarchicalPrior = new GammaPrior(1,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		// TODO Abstract this somewhat
	}

	public void addChildMove(ContinuousPositiveStructAlignMove child) {
		children.add(child);
		//child.addParent(this);
	}
	public double getLogChildDensity(ContinuousPositiveStructAlignMove child) {
		return hierarchicalPrior.logDensity(child.getParam().get());
	}
	public void onlySampleIfAtLeastTwoChildrenNotFixed() {
		onlySampleIfAtLeastTwoChildrenNotFixed = true;
	}
	public void alwaysSample() {
		onlySampleIfAtLeastTwoChildrenNotFixed = false;
	}
	public void setSpikeParams(double[] params) {
		spikeParams = params;
		allowSpikeSelection = true;
		for (ContinuousPositiveStructAlignMove child : children) {
			//System.out.print(child.fixedToParent+" ");
			//double fixProb = spikeParams[0]/(spikeParams[0]+spikeParams[1]);
			//child.fixedToParent = (Utils.generator.nextDouble() < fixProb);
			child.fixedToParent = true;
			//System.out.print("("+child.fixedToParent+") ");
		}
	}
	int countFixedToParent() {
		int n = 0;
		for (ContinuousPositiveStructAlignMove child : children) {
			if (child.fixedToParent) n++;
			//System.out.print(child.fixedToParent+" ");
		}
		//System.out.println();
		return n;
	}
	int countChildren() {
		return children.size();
	}
	@Override
	public double proposal(Object externalState) {		
		if (onlySampleIfAtLeastTwoChildrenNotFixed && (countChildren() - countFixedToParent() < 2)){
			if (children.get(0).parentPriors.get(0) != this) {				
				// Then this must be a nuMove
				// and in this case we keep nu the same,
				// because all (or all but one) of the local sigmas are fixed at 
				// the global.
				autoTune = false;
				proposalWidthControlVariable = 0.1;
				return 0;
			}			
		}
		autoTune = true;
		
		double logProposalDensity = super.proposal(externalState);
		
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {			
				continue;
			}
			if (children.get(i).fixedToParent && children.get(i).parentPriors.get(0) == this) {
				children.get(i).copyState(externalState);
				logProposalDensity -= children.get(i).logPriorDensity(externalState);
			}
			else logProposalDensity -= hierarchicalPrior.logDensity(children.get(i).getParam().get());
		}
		//hierarchicalPrior = new GammaPrior(((StructAlign) owner).nu,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		hierarchicalPrior = new LogNormalPrior(Math.log(((StructAlign) owner).sigma2Hier),((StructAlign) owner).nu);
		//hierarchicalPrior = new GammaPrior(1,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		// TODO Abstract this somewhat
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {
				continue;
			}
			if (children.get(i).fixedToParent && children.get(i).parentPriors.get(0) == this) {
				children.get(i).setParam(param.get());
				logProposalDensity += children.get(i).logPriorDensity(externalState);
			}
			else logProposalDensity += hierarchicalPrior.logDensity(children.get(i).getParam().get());
			//children.get(i).setPlottable();
		}
		return logProposalDensity;
	}
	
	@Override
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {
				continue;
			}
			if (children.get(i).fixedToParent && children.get(i).parentPriors.get(0) == this) {
				children.get(i).restoreState(externalState);				
			}
		}		
		//hierarchicalPrior = new GammaPrior(((StructAlign) owner).nu,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		hierarchicalPrior = new LogNormalPrior(Math.log(((StructAlign) owner).sigma2Hier),((StructAlign) owner).nu);
		//hierarchicalPrior = new GammaPrior(1,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		// TODO Abstract this somewhat
	}
}
