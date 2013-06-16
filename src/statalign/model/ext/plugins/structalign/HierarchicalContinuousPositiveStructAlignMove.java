package statalign.model.ext.plugins.structalign;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Utils;
import statalign.mcmc.GammaPrior;
import statalign.mcmc.McmcMove;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class HierarchicalContinuousPositiveStructAlignMove extends ContinuousPositiveStructAlignMove {

	private List<ContinuousPositiveStructAlignMove> children = new ArrayList<ContinuousPositiveStructAlignMove>();
	public PriorDistribution<Double> hierarchicalPrior;
	
	boolean allowSpikeSelection = false;
	public void disallowSpikeSelection() { allowSpikeSelection = false; }
	public void allowSpikeSelection() { allowSpikeSelection = true; }
	double[] spikeParams;
	
	public HierarchicalContinuousPositiveStructAlignMove (StructAlign s, 
			ParameterInterface p,  PriorDistribution<Double> pr,
			ProposalDistribution<Double> prop, String n) {
		super(s,p,pr,prop,n);
		hierarchicalPrior = new GammaPrior(((StructAlign) owner).nu,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		// TODO Abstract this somewhat
	}

	public void addChildMove(ContinuousPositiveStructAlignMove child) {
		children.add(child);
		child.addParent(this);
	}
	public double getLogChildDensity(ContinuousPositiveStructAlignMove child) {
		return hierarchicalPrior.logDensity(child.getParam().get());
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
		double logProposalDensity = 0;
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {			
				continue;
			}
			if (children.get(i).fixedToParent) {
				children.get(i).copyState(externalState);
				logProposalDensity -= children.get(i).logPriorDensity(externalState);
			}
			else logProposalDensity -= hierarchicalPrior.logDensity(children.get(i).getParam().get());
		}
		logProposalDensity += super.proposal(externalState);
		hierarchicalPrior = new GammaPrior(((StructAlign) owner).nu,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier); 
		// TODO Abstract this somewhat
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {
				continue;
			}
			if (children.get(i).fixedToParent) {
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
			if (children.get(i).fixedToParent) {
				children.get(i).restoreState(externalState);				
			}
		}
		hierarchicalPrior = new GammaPrior(((StructAlign) owner).nu,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		// TODO Abstract this somewhat
	}
}
