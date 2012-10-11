package statalign.model.ext;

import java.io.File;
import java.util.List;

import javax.swing.JComponent;

import statalign.base.Tree;
import statalign.base.Vertex;
import statalign.base.hmm.Hmm;
import statalign.io.DataType;
import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.PluginParameters;

/**
 * Ancestral class for model extension plugins.
 * 
 * @author novak
 */
public abstract class ModelExtension {
	
	public boolean active;
	
	/**
	 * Called during StatAlign startup, after all ModelExtension plugins have been loaded.
	 * 
	 * @param manager reference to manager, save for future use
	 * @param params 
	 */
	public void init(ModelExtManager manager, PluginParameters params) {}
	
	public List<JComponent> getToolBarItems() {
		return null;
	}

	public void setActive(boolean active) {
		this.active = active;
	}
	
	
	public void dataAdded(File file, DataType data) {}

	/**
	 * This should return the log of the model's contribution to the likelihood, it will be added on to
	 * the log-likelihood of the current point in the MCMC state space. Normally will be called in each
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
	 * subalignment have been selected.
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
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {}
	
	public void beforeIndelParamChange(Tree tree, Hmm hmm, int ind) {}
	public void afterIndelParamChange(Tree tree, Hmm hmm, int ind, boolean accepted) {}
	
	public void beforeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {}
	public void afterSubstParamChange(Tree tree, SubstitutionModel model, int ind, boolean accepted) {}
}
