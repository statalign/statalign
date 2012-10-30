package statalign.model.ext;

import java.io.File;
import java.util.List;

import javax.swing.JComponent;

import statalign.base.InputData;
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
	
	boolean active;

	/**
	 * @return <code>true</code> if this plugin is active
	 */
	public boolean isActive() {
		return active;
	}

	/**
	 * Enables or disables this plugin.
	 * @param active <code>true</code> if plugin should be enabled
	 */
	public void setActive(boolean active) {
		this.active = active;
	}

	/**
	 * Override this to return the list of toolbar items to be added. By default returns <code>null</code>.
	 */
	public List<JComponent> getToolBarItems() {
		return null;
	}
	
	public void dataAdded(File file, DataType data) {}

	/**
	 * Called during StatAlign startup, after all ModelExtension plugins have been loaded and command line
	 * has been processed (if present).
	 * 
	 * <p>It is the plugin's responsibility to activate itself using {@link #setActive(boolean)}
	 * when the appropriate command line arguments have been specified.
	 * 
	 * @param manager reference to manager, save for future use
	 * @param params reference to plugin command line parameters or null when StatAlign was run with a GUI
	 */
	public void init(ModelExtManager manager, PluginParameters params) {}
	
	/**
	 * Called during the initialisation of a run if plugin is active. Can be used to process input data.
	 * 
	 * @param input the input data
	 * @exception IllegalArgumentException if the run should be terminated as input data is illegal
	 */
	public void initRun(InputData inputData) throws IllegalArgumentException {}
	
	/**
	 * Called before the start of MCMC sampling, but after the initial tree, alignment etc. have been
	 * generated. Can be used to initialise data structures etc.
	 * @param tree the starting tree
	 */
	public void beforeSampling(Tree tree) {}
	
	public void afterSampling() {}
	
	/**
	 * This should return the log of the model's contribution to the likelihood, it will be added on to
	 * the log-likelihood of the current point in the MCMC state space. Normally will be called in each
	 * MCMC step, when proposing any change.
	 * @param tree current tree
	 * @return log of model extension likelihood, conditional on current tree, alignment and params
	 */
	public abstract double logLikeFactor(Tree tree);
	
	/**
	 * This should return the log of the total prior calculated for the model parameters. It is only used
	 * in parallel mode when proposing swaps between chains. By default returns 0.
	 */
	public double logPrior() {
		return 0;
	}
	
	/**
	 * Returns the weight for choosing a parameter change for this model extension in the MCMC kernel.
	 * By default, returns 0, preventing {@link #proposeParamChange(Tree)} from ever being called.
	 */
	public int getParamChangeWeight() {
		return 0;
	}
	
	/**
	 * Called when this plugin was selected to attempt a model parameter change. The return value
	 * must be the difference between {@link #logLikeFactor(Tree)} after and before
	 * the change.
	 * 
	 * @param tree the current tree
	 * @return the (signed) change in model log-likelihood as a result of the parameter change
	 */
	public double proposeParamChange(Tree tree) {
		return 0;
	}

	/**
	 * Called before an alignment change is proposed, but after the affected subtree has been selected.
	 * May later change to be called after subalignment (window) has also been selected.
	 * @param tree the current tree
	 * @param selectRoot root of the selected subtree
	 */
	public void beforeAlignChange(Tree tree, Vertex selectRoot) {}
	
	/**
	 * Called after an alignment change proposal (accepted or rejected).
	 * @param tree the current tree
	 * @param selectRoot root of the selected subtree
	 * @param accepted true if the change was accepted
	 */
	public void afterAlignChange(Tree tree, Vertex selectRoot, boolean accepted) {}
	
	/**
	 * Called before a topology change is proposed, but after the affected branches are
	 * selected.
	 * @param tree the current tree
	 * @param nephew node that is proposed to be swapped with its "uncle" 
	 */
	public void beforeTreeChange(Tree tree, Vertex nephew) {}
	
	/**
	 * Called after a proposed topology change (accepted or rejected).
	 * @param tree the tree after the change
	 * @param nephew the new nephew (this is the deepest lying node that was affected by the change if there was one)
	 * @param accepted true if the change was accepted
	 */
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {}

	/**
	 * Called before an edge length change is proposed, but after the affected edge is selected.
	 * @param tree the current tree
	 * @param vertex the node whose edge to its parent is selected to be changed
	 */
	public void beforeEdgeLenChange(Tree tree, Vertex vertex) {}
	
	/**
	 * Called after a proposed edge length change (accepted or rejected).
	 * @param tree the current tree
	 * @param vertex the node whose edge to its parent was selected
	 * @param accepted true if the change was accepted
	 */
	public void afterEdgeLenChange(Tree tree, Vertex vertex, boolean accepted) {}
	
	/**
	 * Called before an indel parameter change is proposed, but after the affected parameter is selected.
	 * @param tree the current tree
	 * @param hmm the TKF92 HMM containing its parameters
	 * @param ind the index of the selected parameter
	 */
	public void beforeIndelParamChange(Tree tree, Hmm hmm, int ind) {}
	
	/**
	 * Called after a proposed indel parameter change (accepted or rejected).
	 * @param tree the current tree
	 * @param hmm the TKF92 HMM containing its parameters
	 * @param ind the index of the selected parameter
	 * @param accepted true if the change was accepted
	 */
	public void afterIndelParamChange(Tree tree, Hmm hmm, int ind, boolean accepted) {}

	/**
	 * Called before a substitution parameter change is proposed.
	 * @param tree the current tree
	 * @param model the active substitution model
	 * @param ind the index of the substitution parameter selected to be changed or -1 if unknown 
	 */
	public void beforeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {}
	
	/**
	 * Called after a proposed indel parameter change (accepted or rejected).
	 * @param tree the current tree
	 * @param model the active substitution model
	 * @param ind the index of the substitution parameter selected to be changed or -1 if unknown 
	 * @param accepted true if the change was accepted
	 */
	public void afterSubstParamChange(Tree tree, SubstitutionModel model, int ind, boolean accepted) {}
}
