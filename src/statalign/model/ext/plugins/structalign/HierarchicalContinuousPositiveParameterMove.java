package statalign.model.ext.plugins.structalign;

import java.util.List;

import statalign.base.Tree;
import statalign.model.ext.PriorDistribution;
import statalign.model.ext.plugins.StructAlign;
import statalign.model.ext.plugins.structalign.StructAlignParameterInterface.*;

public class HierarchicalContinuousPositiveParameterMove extends ContinuousPositiveParameterMove {

	private List<ContinuousPositiveParameterMove> children;
	
	public HierarchicalContinuousPositiveParameterMove (StructAlign s, 
			ParameterInterface p, PriorDistribution<Double> pr, 
			String n, double a, double b) {
		super(s,p,pr,n,a,b);
	}
	public HierarchicalContinuousPositiveParameterMove (StructAlign s, 
			ParameterInterface p,  PriorDistribution<Double> pr, 
			String n) {
		super(s,p,pr,n);
	}

	public void addChildMove(ContinuousPositiveParameterMove child) {
		children.add(child);
	}
	public double proposal(Object externalState) {
		double logProposalDensity = super.proposal(externalState);
		// The super method also acquires the Tree
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {
				continue;
			}
			logProposalDensity -= children.get(i).logPriorDensity(externalState);
		}
		owner.sigma2Prior.updateDistribution(owner.nu * owner.sigma2Hier,owner.nu);
		for (int i=0; i<children.size(); i++) {
			if(i == tree.root.index) {
				continue;
			}
			logProposalDensity += children.get(i).logPriorDensity(externalState);
		}
		return logProposalDensity;
	}
	public double logPriorDensity(Object externalState) {
		if (param.get() < minValue) {
			return(Double.NEGATIVE_INFINITY);
		}
		else {
			return prior.logDensity(param.get());
		}
	}
	
	public void updateLikelihood(Object externalState) {
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		else {
			throw new IllegalArgumentException("ContinuousPositiveParameterMove.updateLikelihood must take an argument of type Tree.");
		}
		owner.fullCovar = owner.calcFullCovar(tree);
		owner.curLogLike = owner.calcAllColumnContrib();
	}
	public void restoreState() {
		param.set(oldpar);
		owner.sigma2Prior.updateDistribution(owner.nu * owner.sigma2Hier,owner.nu);
		owner.fullCovar = oldcovar;
		owner.curLogLike = oldll;
	}
}
