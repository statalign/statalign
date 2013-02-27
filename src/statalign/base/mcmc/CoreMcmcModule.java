package statalign.base.mcmc;

import statalign.base.Tree;
import statalign.mcmc.McmcModule;
import statalign.model.ext.ModelExtManager;

public class CoreMcmcModule extends McmcModule {

	ModelExtManager modelExtMan;

	public double totalLogPrior(Tree tree) {
		return(tree.getLogPrior());
	}
	public double logLikeFactor(Tree tree) {
		return(tree.getLogLike());
	}
	@Override
	public boolean isParamChangeAccepted(double logProposalRatio) {
		return false;
	}

}
