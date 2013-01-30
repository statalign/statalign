package statalign.model.ext.plugins.structalign;

import java.util.List;
import java.util.ArrayList;

import statalign.model.ext.McmcMove;
import statalign.model.ext.GammaPrior;
import statalign.model.ext.ParameterInterface;
import statalign.model.ext.PriorDistribution;
import statalign.model.ext.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class HierarchicalContinuousPositiveParameterMove extends ContinuousPositiveParameterMove {

	private List<McmcMove> children = new ArrayList<McmcMove>();
	public PriorDistribution<Double> hierarchicalPrior;
	
	public HierarchicalContinuousPositiveParameterMove (StructAlign s, 
			ParameterInterface p,  PriorDistribution<Double> pr,
			ProposalDistribution<Double> prop, String n) {
		super(s,p,pr,prop,n);
		hierarchicalPrior = new GammaPrior(owner.nu * owner.sigma2Hier,owner.nu);;
		// TODO Abstract this somewhat
	}

	public void addChildMove(ContinuousPositiveParameterMove child) {
		children.add(child);
		child.addParent(this);
	}
	public double getLogChildDensity(ContinuousPositiveParameterMove child) {
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
			logProposalDensity -= children.get(i).logPriorDensity(externalState);
		}
		hierarchicalPrior = new GammaPrior(owner.nu * owner.sigma2Hier,owner.nu); 
		// TODO Abstract this somewhat
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {
				//children.get(i).unsetPlottable();
				continue;
			}
			//children.get(i).setPlottable();
			logProposalDensity += children.get(i).logPriorDensity(externalState);
		}
		return logProposalDensity;
	}
	
	@Override
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		hierarchicalPrior = new GammaPrior(owner.nu * owner.sigma2Hier,owner.nu);
		// TODO Abstract this somewhat
	}
}
