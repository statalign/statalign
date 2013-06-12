package statalign.model.ext.plugins.structalign;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.mcmc.ContinuousPositiveParameterMove;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class ContinuousPositiveStructAlignMove extends ContinuousPositiveParameterMove {

	StructAlign structAlign;
	double[][] oldcovar;
	boolean oldFixedToParent;
	
	List<HierarchicalContinuousPositiveStructAlignMove> parentPriors = null;
	
	/** 
	 * <code>true</code> if this parameter is currently constrained to be
	 * equal to its parent (where uniquely defined). This is used in cases
	 * such as spike and slab priors.  
	 */
	boolean fixedToParent = false;
	int nFixedToParent = 0;

	public StructAlignMoveParams moveParams = new StructAlignMoveParams();
	
	public void addParent(HierarchicalContinuousPositiveStructAlignMove p) {
		if (parentPriors == null) {
			parentPriors = new ArrayList<HierarchicalContinuousPositiveStructAlignMove>();
		}
		parentPriors.add(p);
	}
	
	public ContinuousPositiveStructAlignMove (StructAlign s, 
			ParameterInterface p, PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		super(s,p,pr,prop,n);
		structAlign = s;
	}
	public void copyState(Object externalState) {
		super.copyState(externalState);
		oldcovar = structAlign.fullCovar;
		oldFixedToParent = fixedToParent;
	}
	
	@Override
	public double proposal(Object externalState) {		
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("ContinuousPositiveStructAlignMove.proposal must take an argument of type Tree.");
		}
		proposalDistribution.updateProposal(proposalWidthControlVariable,param.get());
		
		double logProposalDensity = 0;
		
		if (parentPriors != null) {
			for (HierarchicalContinuousPositiveStructAlignMove parent : parentPriors) {
				if (parent.allowSpike) { 			
					double[] spikeParams = parent.spikeParams;
					int m = parent.countFixedToParent();
					int n = parent.countChildren();
					double fixProb = spikeParams[0]/(spikeParams[0]+spikeParams[1]);
					fixedToParent = (Utils.generator.nextDouble()<fixProb);
					if (fixedToParent) {
						param.set(parent.param.get());					 
						if (!oldFixedToParent) { // z_k = 0 to begin with
							logProposalDensity += Math.log((1-fixProb)/fixProb);
							logProposalDensity += Math.log((m+1)/(n-m) 
											 	* (spikeParams[0]+m) 
											 	/ (spikeParams[1]+n-m-1));
						}
					}
					else {
						if (oldFixedToParent) { // z_k = 1 to begin with
							logProposalDensity += Math.log(fixProb/(1-fixProb));
							logProposalDensity += Math.log((n-m+1)/m 
								 				* (spikeParams[1]+n-m) 
								 				/ (spikeParams[0]+m-1));
						}
					}
					// Assume only one parent allows spikes, and break once we've
					// found it.
					break;
				}
			}
		}
		moveParams.setFixedToParent(fixedToParent);
		if (fixedToParent) nFixedToParent++;
		else param.set(proposalDistribution.sample());
		
		if (param.get() < minValue || param.get() > maxValue) {
			return(Double.NEGATIVE_INFINITY);
		}

		/** - log p(new | old) */
		if (!fixedToParent) logProposalDensity = -proposalDistribution.logDensity(param.get());
		
		proposalDistribution.updateProposal(proposalWidthControlVariable,param.get());
		
		/** + log p(old | new) */
		if (!oldFixedToParent) logProposalDensity += proposalDistribution.logDensity(oldpar);
		
		return logProposalDensity;
	}
	@Override
	public double logPriorDensity(Object externalState) {
		if (param.get() < minValue) {
			return(Double.NEGATIVE_INFINITY);
		}
		else {
			// return 0;
			double logDensity = prior.logDensityUnnormalised(param.get());
			// Since we're only using this in ratios, there's no
			// need to compute the normalising constant, which is good, 
			// because some priors may be improper.
			// NB be careful with this though -- an improper prior should
			// only be used if the posterior can be shown to be proper.
			if (parentPriors != null) {
				for (HierarchicalContinuousPositiveStructAlignMove parent : parentPriors) {
					logDensity += parent.getLogChildDensity(this);
					// The normalising constant of this density will depend
					// on the parent, so we may need the normalised density here.
				}
			}
			return logDensity;
		}
	}
	public void updateLikelihood(Object externalState) {
		if (param.get() > minValue) {
			structAlign.fullCovar = structAlign.calcFullCovar(tree);
			owner.setLogLike( structAlign.calcAllColumnContrib() );
		}
	}
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		structAlign.fullCovar = oldcovar;
		fixedToParent = oldFixedToParent;
		moveParams.setFixedToParent(fixedToParent);
	}
	public void setParam(double x) {
		param.set(x);
	}
	
}
