package statalign.model.ext.plugins.structalign;

import java.util.List;
import java.util.ArrayList;
import statalign.base.mcmc.ContinuousPositiveParameterMove;
import statalign.base.mcmc.ParameterInterface;
import statalign.base.mcmc.PriorDistribution;
import statalign.base.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class ContinuousPositiveStructAlignMove extends ContinuousPositiveParameterMove {

	StructAlign structAlign;
	double[][] oldcovar;
	
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
		oldcovar = structAlign.fullCovar;
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
	}
}
