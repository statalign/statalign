package statalign.model.ext.plugins.structalign;

import java.util.List;
import java.util.ArrayList;
import statalign.model.ext.ParameterInterface;
import statalign.model.ext.PriorDistribution;
import statalign.model.ext.ProposalDistribution;
import statalign.base.Tree;
import statalign.model.ext.plugins.StructAlign;
import statalign.utils.GammaDistribution;
import statalign.model.ext.plugins.structalign.StructAlignParameterInterface.*;

public class ContinuousPositiveParameterMove extends StructAlignMcmcMove {

	Tree tree;
	List<HierarchicalContinuousPositiveParameterMove> parentPriors = null;
	
	public void addParent(HierarchicalContinuousPositiveParameterMove p) {
		if (parentPriors == null) {
			parentPriors = new ArrayList<HierarchicalContinuousPositiveParameterMove>();
		}
		parentPriors.add(p);
	}

	private ProposalDistribution<Double> proposalDistribution;

	double oldpar;
	double[][] oldcovar;
	double oldll;
	
	protected double minValue = 0.0;
	// If the proposal takes the parameter below this value, 
	// the move is rejected.
	
	public ContinuousPositiveParameterMove (StructAlign s, 
			ParameterInterface p, PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		owner = s;
		param = p;
		prior = pr;
		name = n;
		proposalDistribution = prop;
		autoTune = true;
	}
	public void setMinValue(double x) {
		minValue = x;
	}
	public void copyState(Object externalState) {
		oldpar = param.get();
		oldcovar = owner.fullCovar;
		oldll = owner.curLogLike;
	}

	public double proposal(Object externalState) {

		proposalDistribution.updateProposal(proposalWidthControlVariable,param.get());
		param.set(proposalDistribution.sample());
		
		if (param.get() < minValue) {
			return(Double.NEGATIVE_INFINITY);
		}

		/** - log p(new | old) */
		double logProposalDensity = -proposalDistribution.logDensity(param.get());
		
		proposalDistribution.updateProposal(proposalWidthControlVariable,param.get());
		
		/** + log p(old | new) */
		logProposalDensity += proposalDistribution.logDensity(oldpar);
		
		return logProposalDensity;
	}
	
	public double logPriorDensity(Object externalState) {
		if (param.get() < minValue) {
			return(Double.NEGATIVE_INFINITY);
		}
		else {
			// return 0;
			double logDensity = prior.logDensityUnnormalised(param.get());
			// Since we're only using this in ratios, there's no
			// need to compute the normalising constant.
			if (parentPriors != null) {
				for (HierarchicalContinuousPositiveParameterMove parent : parentPriors) {
					logDensity += parent.getLogChildDensity(this);
				}
			}
			return logDensity;
		}
	}
	public void updateLikelihood(Object externalState) {
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		else {
			throw new IllegalArgumentException("ContinuousPositiveParameterMove.updateLikelihood must take an argument of type Tree.");
		}
		if (param.get() > minValue) {
			((StructAlign) owner).fullCovar = owner.calcFullCovar(tree);
			owner.curLogLike = owner.calcAllColumnContrib();
		}
	}
	public void restoreState(Object externalState) {
		param.set(oldpar);
		owner.fullCovar = oldcovar;
		owner.curLogLike = oldll;
	}
}
