package statalign.base.mcmc;

import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.mcmc.McmcModule;
import statalign.model.ext.ModelExtManager;

public class CoreMcmcModule extends McmcModule {

	private ModelExtManager modelExtMan;
	public ModelExtManager getModelExtMan() {
		return modelExtMan;
	}
	
	public CoreMcmcModule (Mcmc mc, ModelExtManager m) {
		setMcmc(mc);
		modelExtMan = m;
	}

	public double totalLogPrior(Tree tree) {
		return(tree.getLogPrior());
	}
	public double logLikeFactor(Tree tree) {
		return(tree.getLogLike());
	}
//	@Override
//	public boolean isParamChangeAccepted(double logProposalRatio) {
//		return mcmc.modExtParamChangeCallback(logProposalRatio);
//	}
//	@Override
//	public boolean isParamChangeAccepted(double logProposalRatio) {
//		return mcmc.isParamChangeAccepted(logProposalRatio);
//	}

}
