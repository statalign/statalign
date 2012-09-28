package statalign.model.ext;

import statalign.base.Tree;
import statalign.base.Vertex;
import statalign.base.hmm.Hmm;
import statalign.model.subst.SubstitutionModel;

/**
 * Ancestral class for model extension plugins.
 * 
 * @author novak
 */
public abstract class ModelExtension {

	/**
	 * This should return the log of the model's contribution to the likelihood, it will be added on to
	 * the likelihood of the current point in the MCMC state space. Normally will be called in each
	 * MCMC step, when proposing any change.
	 * @param tree current tree
	 * @return log of model extension likelihood, conditional on current tree, alignment and params
	 */
	public abstract double logLikeFactor(Tree tree);
	
	/**
	 * Returns the weight for choosing a parameter change for this model extension in the MCMC kernel.
	 * By default, returns 0, preventing {@link #proposeParamChange(Tree, double)} from ever being called.
	 */
	public int getParamChangeWeight() {
		return 0;
	}
	
	/**
	 * Called when
	 * 
	 * @param tree
	 * @param logLike
	 * @return
	 */
	public boolean proposeParamChange(Tree tree, double logLike) {
		return false;
	}

	/**
	 * Called before an alignment change is proposed, but after the affected subtree and
	 * subalignment has been selected.
	 * @param tree the main Tree object
	 * @param selectRoot root of the selected subtree
	 */
	public void beforeAlignChange(Tree tree, Vertex selectRoot) {}
	/**
	 * Called after an alignment change proposal (accepted or rejected).
	 * @param tree the main Tree object
	 * @param selectRoot root of the selected subtree
	 * @param accepted true if the change was accepted
	 */
	public void afterAlignChange(Tree tree, Vertex selectRoot, boolean accepted) {}
	
	/**
	 * Called before a topology change is proposed, but after the affected branches are
	 * selected.
	 * @param tree the main Tree object
	 * @param nephew node that is proposed to be swapped with its uncle 
	 */
	public void beforeTreeChange(Tree tree, Vertex nephew) {}
	/**
	 * Called after a proposed topology change (accepted or rejected).
	 * @param tree
	 * @param nephew
	 */
	public void afterTreeChange(Tree tree, Vertex nephew) {}
	
	public void beforeIndelParamChange(Tree tree, Hmm hmm, int ind) {}
	public void afterIndelParamChange(Tree tree, Hmm hmm, int ind) {}
	
	public void beforeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {}
	public void afterSubstParamChange(Tree tree, SubstitutionModel model, int ind) {}
}
