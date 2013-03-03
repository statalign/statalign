package statalign.model.ext.plugins.structalign;

import java.util.List;
import java.util.ArrayList;

import statalign.mcmc.GammaPrior;
import statalign.mcmc.McmcMove;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class HierarchicalContinuousPositiveStructAlignMove extends ContinuousPositiveStructAlignMove {

	private List<McmcMove> children = new ArrayList<McmcMove>();
	public PriorDistribution<Double> hierarchicalPrior;
	
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
	
	@Override
	public double proposal(Object externalState) {
		double logProposalDensity = super.proposal(externalState);
		// The super method also acquires the Tree
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {
				continue;
			}
			logProposalDensity -= hierarchicalPrior.logDensity(children.get(i).getParam().get());
		}
		hierarchicalPrior = new GammaPrior(((StructAlign) owner).nu,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier); 
		// TODO Abstract this somewhat
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {
				continue;
			}
			//children.get(i).setPlottable();
			logProposalDensity += hierarchicalPrior.logDensity(children.get(i).getParam().get());
		}
		return logProposalDensity;
	}
	
	@Override
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		hierarchicalPrior = new GammaPrior(((StructAlign) owner).nu,((StructAlign) owner).nu / ((StructAlign) owner).sigma2Hier);
		// TODO Abstract this somewhat
	}
}
