package statalign.base.mcmc;

import java.util.List;
import java.util.ArrayList;

import statalign.base.Tree;
import statalign.base.Utils;
import statalign.mcmc.ContinuousPositiveParameterMove;
import statalign.mcmc.McmcModule;
import statalign.mcmc.McmcMove;
import statalign.mcmc.ParameterInterface;
import statalign.mcmc.PriorDistribution;
import statalign.mcmc.ProposalDistribution;
import statalign.model.ext.plugins.StructAlign;

public class SubstMove extends McmcMove {
	
	Tree tree = null;
	
	public SubstMove (McmcModule m, String n) {
		param = null; // SubstitutionModel handles the parameter access itself
		owner = m;
		prior = null;
		name = n;
		proposalWidthControlVariable = 0.1;		
		maxProposalWidthControlVariable = 10.0;
		spanMultiplier = 0.9;
		//autoTune = false;
	}

	public void copyState(Object externalState) {
		if (externalState instanceof Tree) {
			if (tree == null) {
				tree = (Tree) externalState;
			}
		}
		else {
			throw new IllegalArgumentException("SubstMove.move must take an argument of type Tree.");
		}
		if (tree.substitutionModel.params.length == 0) {
			return;
		}
		((CoreMcmcModule) owner).getModelExtMan().beforeSubstParamChange(tree,tree.substitutionModel, -1);
		// SubstitutionModel does the rest internally
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
		tree.substitutionModel.updateSpan(proposalWidthControlVariable);
		double logProposalDensity = -tree.substitutionModel.sampleParameter();
		// For some reason all the SubstitutionModel sampleParameter() functions
		// return the inverse of the Hastings ratio, hence the minus sign above.
		return logProposalDensity;
	}
	
	public double logPriorDensity(Object externalState) {
		return Math.log(tree.substitutionModel.getPrior());
	}
	public void updateLikelihood(Object externalState) {
		for (int i = 0; i < tree.vertex.length; i++) {
			tree.vertex[i].updateTransitionMatrix();
		}
		tree.root.calcFelsenRecursively();
		owner.setLogLike(((CoreMcmcModule) owner).getModelExtMan().logLikeSubstParamChange(tree, tree.substitutionModel, -1));
	}
		
	@Override
	public void afterMove(Object externalState) {
		((CoreMcmcModule) owner).getModelExtMan().afterSubstParamChange(tree, tree.substitutionModel, -1, lastMoveAccepted);
	}
	
	@Override
	public void restoreState(Object externalState) {
		tree.substitutionModel.restoreParameter();
		for (int i = 0; i < tree.vertex.length; i++) {
			tree.vertex[i].updateTransitionMatrix();
		}
		tree.root.calcFelsenRecursively();
	}


}
