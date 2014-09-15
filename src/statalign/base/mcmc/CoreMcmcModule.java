package statalign.base.mcmc;

import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.mcmc.McmcModule;
import statalign.model.ext.ModelExtManager;

/**
 * curLogLike for this McmcModule includes the
 * contributions from the ModelExtensions.
 * 
 * The core McmcMoves need to call functions in modelExtMan
 * before and after execution, and also in order to update 
 * the likelihood contribution from all the ModelExtensions.   
 * 
 * @author herman
 * 
 */
public class CoreMcmcModule extends McmcModule {
	
	private ModelExtManager modelExtMan;
	
	public ModelExtManager getModelExtMan() {
		return modelExtMan;
	}
	
	public CoreMcmcModule (Mcmc mc, ModelExtManager m) {
		super("coreModel");
		setMcmc(mc);
		modelExtMan = m;
	}

	public double totalLogPrior(Tree tree) {		
		return(tree.getLogPrior());
	}
	public double logLikeFactor(Tree tree) {
		return(tree.getLogLike());
	}
	@Override
	public boolean proposeParamChange(Tree tree) {
		// Each move will update the curLogLike variable
		// for `this'
		boolean accepted = false;
		int[] weights = new int[2];
		weights[0] = getParamChangeWeight();
		weights[1] = modelExtMan.getParamChangeWeight();
		int choice = Utils.weightedChoose(weights);
		if (choice == 0) { // Core move selected
			//System.out.println("Core move.");
			accepted = super.proposeParamChange(tree);
		}
		else { // Model extension move selected
			//System.out.println("Model extension move.");
			modelExtMan.beforeModExtParamChange(tree);
			accepted = modelExtMan.proposeParamChange(tree);
			modelExtMan.afterModExtParamChange(tree, accepted);
		}
		return accepted;
	}
	
	@Override
	public void beforeSampling(Tree tree) {
		super.beforeSampling(tree);
		modelExtMan.beforeSampling(tree);
	}
} 
