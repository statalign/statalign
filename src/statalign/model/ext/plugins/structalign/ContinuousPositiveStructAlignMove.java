package statalign.model.ext.plugins.structalign;

import java.util.HashMap;
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
	//private HashMap<Integer, MultiNormCholesky> oldMultiNorms; // TODO handle this inside StructAlign
	//private HashMap<Integer, MultiNormCholesky> oldMultiNormsLocal; // TODO handle this inside StructAlign	
	boolean oldFixedToParent;
	
	List<HierarchicalContinuousPositiveStructAlignMove> parentPriors = null;
	
	/** 
	 * <code>true</code> if this parameter is currently constrained to be
	 * equal to its parent (where uniquely defined). This is used in cases
	 * such as spike and slab priors.  
	 */
	boolean fixedToParent = false;

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

		oldcovar = structAlign.fullCovar; // TODO handle in more abstract fashion
		oldll = structAlign.curLogLike; // TODO handle in more abstract fashion
		
		//oldMultiNorms = structAlign.multiNorms; // TODO handle in more abstract fashion
		structAlign.beforeContinuousParamChange(tree);
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
		boolean allowKeepingOfParameter = false;
		if (parentPriors != null) {
			for (HierarchicalContinuousPositiveStructAlignMove parent : parentPriors) {
				if (parent.allowSpikeSelection) { 			
					double[] spikeParams = parent.spikeParams;
					double m = (double) parent.countFixedToParent();
					double n = (double) parent.countChildren();
					allowKeepingOfParameter = owner.isBurnin() & (n-m < 2); 
					// In this case the hierarchical variance could become very small, 
					// so we need to allow switching of z_k to zero without resampling
					// its parameter, during the burnin (where no adjustment to
					// the M-H ratio is strictly necessary)
					
					double fixProb = spikeParams[0]/(spikeParams[0]+spikeParams[1]);
					fixedToParent = (Utils.generator.nextDouble()<fixProb);
					if (fixedToParent) {
						param.set(parent.param.get());					 
						if (!oldFixedToParent) { // z_k = 0 to begin with
							// Ratio of proposal probabilities for z_k
							logProposalDensity += Math.log((1-fixProb)/fixProb);
							// Then also include ratio of priors for the vector z
							logProposalDensity += Math.log((n-m)/(m+1) 
											 	* (spikeParams[0]+m) 
											 	/ (spikeParams[1]+n-m-1));
						}
					}
					else {
						if (oldFixedToParent) { // z_k = 1 to begin with
							// Ratio of proposal probabilities for z_k
							logProposalDensity += Math.log(fixProb/(1-fixProb));	
							// Then also include ratio of priors for the vector z
							logProposalDensity += Math.log(m/(n-m+1) 
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
		if (!allowKeepingOfParameter || (Utils.generator.nextDouble() < 0.8)) {
			// If allowKeepingOfParameter, we still resample the parameter
			// with 0.8 probability, to cover the cases where the non-spike
			// density is lower at the mode. In the cases where it's higher,
			// and the variance is very low (less likely), then we can 
			// keep the current parameter and just allow switching to the higher 
			// density.
			if (!fixedToParent) param.set(proposalDistribution.sample());
				
			if (param.get() < minValue || param.get() > maxValue) {
				return(Double.NEGATIVE_INFINITY);
			}
		
			/** - log p(new | old) */
			if (!fixedToParent) logProposalDensity -= proposalDistribution.logDensity(param.get());
			
			proposalDistribution.updateProposal(proposalWidthControlVariable,param.get());
			
			/** + log p(old | new) */
			if (!oldFixedToParent) logProposalDensity += proposalDistribution.logDensity(oldpar);				
		}
		
		return logProposalDensity;
	}
	@Override
	public double logPriorDensity(Object externalState) {
		if (param.get() < minValue) {
			return(Double.NEGATIVE_INFINITY);
		}
		else {
			// Contribution from independent prior on this variable
			double logDensity = prior.logDensityUnnormalised(param.get());
			// Since we're only using this in ratios, there's no
			// need to compute the normalising constant, which is good, 
			// because some priors may be improper.
			// NB be careful with this though -- an improper prior should
			// only be used if the posterior can be shown to be proper.

			// Then add the contribution from the hierarchical prior
			if (!fixedToParent && parentPriors != null) {
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
		if (!(fixedToParent && oldFixedToParent) && param.get() > minValue) {
			owner.setLogLike( structAlign.logLikeContinuousParamChange(tree) );			
		}
	}
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		if (!structAlign.globalSigmaSpike || !(fixedToParent && oldFixedToParent)) {
			structAlign.fullCovar = oldcovar; // TODO handle in more abstract fashion
			structAlign.curLogLike = oldll; // TODO handle in more abstract fashion
			//structAlign.multiNorms = oldMultiNorms; // TODO handle in more abstract fashion
			
			//System.out.println(lastMoveAccepted);
			structAlign.afterContinuousParamChange(tree, false);
			// Can't refer to lastMoveAccepted here, because the restoreState
			// function may be called from an McmcCombinationMove.
			
			if (structAlign.globalSigmaSpike) {
				fixedToParent = oldFixedToParent;
				moveParams.setFixedToParent(fixedToParent);
			}
		}
	}
	public void setParam(double x) {
		param.set(x);
	}
	public void afterMove(Object externalState) {
		super.afterMove(externalState);		
		autoTune = !fixedToParent;
		if (fixedToParent) proposalWidthControlVariable = parentPriors.get(0).proposalWidthControlVariable;
	}
}

