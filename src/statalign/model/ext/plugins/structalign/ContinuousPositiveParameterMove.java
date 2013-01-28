package statalign.model.ext.plugins.structalign;

import statalign.model.ext.McmcMove;
import statalign.base.Tree;
import statalign.model.ext.plugins.StructAlign;
import statalign.utils.GammaDistribution;
import statalign.model.ext.plugins.structalign.StructAlignParameterInterface.*;

public class ContinuousPositiveParameterMove extends McmcMove {

	StructAlign owner;
	Tree tree;

	double oldpar;
	double[][] oldcovar;
	double oldll;
	ParameterInterface param;
	
	public ContinuousPositiveParameterMove (StructAlign s, ParameterInterface p) {
		owner = s;
		param = p;
	}
	public void copyState() {
		oldpar = param.get();
		oldcovar = owner.fullCovar;
		oldll = owner.curLogLike;
	}

	public double proposal(Object externalState) {
		
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		GammaDistribution proposal;
		GammaDistribution reverse;
		
		double newpar = 1.0d/proposalWidthControlVariable;
		proposal = new GammaDistribution(newpar + 0.001, oldpar / newpar + 0.001);
		//proposal = new GammaDistribution((0.001+oldpar)/sigma2P, sigma2P);
		param.set(proposal.sample());
		reverse = new GammaDistribution(newpar + 0.001, param.get() / newpar + 0.001);
		
		return(Math.log(reverse.density(oldpar)) 
				- Math.log(proposal.density(param.get())));
	}
	
	public void updateLikelihood(Object externalState) {
		if (externalState instanceof Tree) {
			tree = (Tree) externalState;
		}
		owner.fullCovar = owner.calcFullCovar(tree);
		owner.curLogLike = owner.calcAllColumnContrib();
	}
	public void restoreState() {
		param.set(oldpar);
		owner.fullCovar = oldcovar;
		owner.curLogLike = oldll;
	}
}
