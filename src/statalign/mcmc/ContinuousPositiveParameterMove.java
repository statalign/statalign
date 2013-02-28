package statalign.mcmc;

import java.util.List;
import java.util.ArrayList;
import statalign.base.Tree;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;
import statalign.utils.GammaDistribution;
import statalign.model.ext.plugins.structalign.StructAlignParameterInterface.*;

public abstract class ContinuousPositiveParameterMove extends McmcMove {

	protected Tree tree = null;

	private ProposalDistribution<Double> proposalDistribution;

	protected double oldpar;
	protected double oldll;
	
	protected double minValue = 0.0;
	protected double maxValue = Double.POSITIVE_INFINITY;
	// If the proposal takes the parameter below/above these values, 
	// the move is rejected.
	
	public ContinuousPositiveParameterMove (McmcModule m,
			ParameterInterface p, 
			PriorDistribution<Double> pr, 
			ProposalDistribution<Double> prop, String n) {
		owner = m;
		param = p;
		prior = pr;
		name = n;
		proposalDistribution = prop;
		autoTune = true;
	}
	public void setMinValue(double x) {
		minValue = x;
	}
	public void setMaxValue(double x) {
		maxValue = x;
	}
	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("ContinuousPositiveParameterMove.copyState must take an argument of type Tree.");
		}
		oldpar = param.get();
		oldll = owner.getLogLike();
	}

	public double proposal(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("ContinuousPositiveParameterMove.proposal must take an argument of type Tree.");
		}
		proposalDistribution.updateProposal(proposalWidthControlVariable,param.get());
		param.set(proposalDistribution.sample());
		
		if (param.get() < minValue || param.get() > maxValue) {
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
		if (param.get() < minValue || param.get() > maxValue) {
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
			return logDensity;
		}
	}
	public abstract void updateLikelihood(Object externalState);
	public void restoreState(Object externalState) {
		param.set(oldpar);
		owner.setLogLike(oldll);
	}
}
