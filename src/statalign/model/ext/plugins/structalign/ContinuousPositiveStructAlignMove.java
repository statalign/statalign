package statalign.model.ext.plugins.structalign;

import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

import statalign.mcmc.ContinuousPositiveParameterMove;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class ContinuousPositiveStructAlignMove extends ContinuousPositiveParameterMove {

	StructAlign structAlign;
	double[][] oldcovar; // TODO handle this inside StructAlign
	private HashMap<Integer, MultiNormCholesky> oldMultiNorms; // TODO handle this inside StructAlign	

	
	List<HierarchicalContinuousPositiveStructAlignMove> parentPriors = null;

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
		oldMultiNorms = structAlign.multiNorms; // TODO handle in more abstract fashion
		structAlign.beforeContinuousParamChange(tree);
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
//			structAlign.fullCovar = structAlign.calcFullCovar(tree);
//			owner.setLogLike( structAlign.calcAllColumnContrib() );
			owner.setLogLike( structAlign.logLikeContinuousParamChange(tree) );			
		}
	}
	public void restoreState(Object externalState) {
		super.restoreState(externalState);
		structAlign.fullCovar = oldcovar; // TODO handle in more abstract fashion
		structAlign.curLogLike = oldll; // TODO handle in more abstract fashion
		structAlign.multiNorms = oldMultiNorms; // TODO handle in more abstract fashion
		structAlign.afterContinuousParamChange(tree, lastMoveAccepted);
	}
}
