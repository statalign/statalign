package statalign.model.ext.plugins.structalign;

import statalign.model.ext.GammaPrior;
import statalign.model.ext.McmcMove;
import statalign.model.ext.PriorDistribution;
import statalign.base.Tree;
import statalign.model.ext.plugins.StructAlign;
import statalign.utils.GammaDistribution;
import statalign.model.ext.plugins.structalign.StructAlignParameterInterface.*;

public class ContinuousPositiveParameterMove extends McmcMove {

	StructAlign owner;
	Tree tree;
	
	/**
	 * Controls the shape of the proposal distribution, relative to
	 * the current parameter value.
	 */
	private double proposalShape;
	private double proposalRate;

	double oldpar;
	double[][] oldcovar;
	double oldll;
	ParameterInterface param;
	
	protected double minValue = 0.0;
	// If the proposal takes the parameter below this value, 
	// the move is rejected.
	
	public ContinuousPositiveParameterMove (StructAlign s, 
			ParameterInterface p, GammaPrior pr, 
			String n, double a, double b) {
		owner = s;
		param = p;
		prior = pr;
		name = n;
		proposalShape = a;
		proposalRate = b;
		autoTune = true;
	}
	public ContinuousPositiveParameterMove (StructAlign s, 
			ParameterInterface p, GammaPrior pr, String n) {
		owner = s;
		param = p;
		prior = pr;
		name = n;
		proposalShape = 0.001;
		proposalRate = 0.001;
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
		
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		else {
			throw new IllegalArgumentException("ContinuousPositiveParameterMove.proposal must take an argument of type Tree.");
		}
		
		if (param instanceof Sigma2Interface) {
			if (((Sigma2Interface) param).getIndex() == tree.root.index) {
				return(Double.NEGATIVE_INFINITY);
			}
		}
		GammaDistribution proposal;
		GammaDistribution reverse;
		
		double conc = 1.0d/proposalWidthControlVariable;
		proposal = new GammaDistribution(conc + proposalShape, 
				oldpar / (conc + proposalRate) );
		param.set(proposal.sample());
		reverse = new GammaDistribution(conc + proposalShape, 
				param.get() / (conc + proposalRate) );
		
		if (param.get() < minValue) {
			return(Double.NEGATIVE_INFINITY);
		}
		else {
			return(Math.log(reverse.density(oldpar)) 
				- Math.log(proposal.density(param.get())));
		}
	}
	public double logPriorDensity(Object externalState) {
		return 0;
	}
	public void updateLikelihood(Object externalState) {
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		else {
			throw new IllegalArgumentException("ContinuousPositiveParameterMove.updateLikelihood must take an argument of type Tree.");
		}
		((StructAlign) owner).fullCovar = owner.calcFullCovar(tree);
		owner.curLogLike = owner.calcAllColumnContrib();
	}
	public void restoreState(Object externalState) {
		param.set(oldpar);
		owner.fullCovar = oldcovar;
		owner.curLogLike = oldll;
	}
}
