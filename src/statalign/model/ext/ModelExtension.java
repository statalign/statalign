package statalign.model.ext;

import java.io.File;
import java.util.List;

import javax.swing.JComponent;

import statalign.base.InputData;
import statalign.base.Mcmc;
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
	
	/** Private access to model extension manager */
	private ModelExtManager manager;
	
	protected boolean active;
	
	public int[] proposalCounts;
	public int[] acceptanceCounts;
	public double[] proposalWidthControlVariables;
	
	public void setManager(ModelExtManager manager) {
		this.manager = manager;
	}

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
	 * Returns the (unmodifiable) list of recognised {@link ModelExtension} plugins.
	 */
	public final List<ModelExtension> getPluginList() {
		return manager.getPluginList();
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
	 * @param params reference to plugin command line parameters or null when StatAlign was run with a GUI
	 */
	public void init(PluginParameters params) {}
	
	/**
	 * Called during the initialisation of a run if plugin is active. Override to process input data.
	 * 
	 * @param input the input data
	 * @exception IllegalArgumentException if the run should be terminated as input data is illegal
	 */
	public void initRun(InputData inputData) throws IllegalArgumentException {}
	
	/**
	 * Called before the start of MCMC sampling, but after the initial tree, alignment etc. have been
	 * generated. Override to initialise data structures etc.
	 * @param tree the starting tree
	 */
	public void beforeSampling(Tree tree) {}
	
	public void afterSampling() {}
	
	/**
	 * This should return the log of the model's contribution to the likelihood, it will be added on to
	 * the log-likelihood of the current point in the MCMC state space. Normally it will be called once at the
	 * initialisation of the MCMC process and from then on once in each MCMC step, when proposing any change.
	 * In debug mode, will be called more often (including after proposed changes) to ensure consistency.
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
	 * Called when a model extension plugin is to propose a parameter change. Plugins must check whether
	 * they are selected for the parameter change and propose a change if necessary. Will be called
	 * for the selected plugin first. {@link #logLikeFactor(Tree)} will be called afterwards to get
	 * the likelihood contribution after the change, and the move will be accepted or rejected by the
	 * framework. All plugins will be notified about its outcome by
	 * {@link #afterModExtParamChange(Tree, ModelExtension, boolean)}.
	 * @param tree the current tree
	 * @param ext the model extension plugin that is selected to propose a parameter change
	 */
	public void beforeModExtParamChange(Tree tree, ModelExtension ext) {}
	
	/**
	 * Called when this plugin was selected to attempt a model parameter change. Plugin should call
	 * {@link #isParamChangeAccepted(double)} with the log of a Metropolis-Hastings-like ratio (see the
	 * method for details) to find out whether the change was accepted and then act
	 * accordingly, ie. keep or restore the state.
	 * 
	 * <p>The call to the above method will result in a call to {@link #logLikeModExtParamChange(Tree, ModelExtension)}
	 * with ModelExtension set to <code>this</code> and the plugin should return the log-likelihood factor
	 * calculated after the proposed change. 
	 * 
	 * <p>All plugins will later be notified about the acceptance/rejection of the change through
	 * {@link #afterModExtParamChange(Tree, ModelExtension, boolean)}.
	 * 
	 * @param tree the current tree
	 */
	public void proposeParamChange(Tree tree) {}
	
	/**
	 * Should be called from {@link #proposeParamChange(Tree)} to find out whether a proposed parameter
	 * change was accepted.
	 * 
	 * <p>Parameter <code>logProbRatio</code> must be the log of P(x|x')/P(x'|x) * Pr(x')/Pr(x) where x is the
	 * old value of the model parameter, x' is new value, P(x'|x) is the probability of the proposed change
	 * (proposal probability), P(x|x') is the backproposal probability, Pr(x') is the prior probability
	 * of the new parameter value (new prior), Pr(x) is the old prior. The remaining factor
	 * Pi(new state)/Pi(old state) of the Metropolis-Hastings ratio will be calculated by calls to
	 * {@link #logLikeModExtParamChange(Tree, ModelExtension)}.
	 * 
	 * @param logProbRatio MH-like ratio as explained above
	 * @return <code>true</code> if the change was accepted
	 */
	public final boolean isParamChangeAccepted(double logProbRatio) {
		return manager.modExtParamChangeCallback(logProbRatio);
	}
	
	public double logLikeModExtParamChange(Tree tree, ModelExtension ext) {
		return logLikeFactor(tree);
	}
	
	/**
	 * Called after a proposed model extension parameter change (accepted or rejected).
	 * @param tree the current tree
	 * @param ext the model extension plugin that was selected to propose a parameter change
	 * @param accepted <code>true</code> if the change was accepted
	 */
	public void afterModExtParamChange(Tree tree, ModelExtension ext, boolean accepted) {}
	
	/**
	 * Called before an alignment change is proposed, but after the affected subtree has been selected.
	 * May change later to be called after subalignment (window) has also been selected.
	 * @param tree the current tree
	 * @param selectRoot root of the selected subtree
	 */
	public void beforeAlignChange(Tree tree, Vertex selectRoot) {}
	
	/**
	 * 
	 * @param tree
	 * @param selectRoot
	 * @return
	 */
	public double logLikeAlignChange(Tree tree, Vertex selectRoot) {
		return logLikeFactor(tree);
	}
	
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
	
	public double logLikeTreeChange(Tree tree, Vertex nephew) {
		return logLikeFactor(tree);
	}
	
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
	
	public double logLikeEdgeLenChange(Tree tree, Vertex vertex) {
		return logLikeFactor(tree);
	}
	
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
	
	public double logLikeIndelParamChange(Tree tree, Hmm hmm, int ind) {
		return logLikeFactor(tree);
	}
	
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
	
	public double logLikeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {
		return logLikeFactor(tree);
	}
	
	/**
	 * Called after a proposed indel parameter change (accepted or rejected).
	 * @param tree the current tree
	 * @param model the active substitution model
	 * @param ind the index of the substitution parameter selected to be changed or -1 if unknown 
	 * @param accepted true if the change was accepted
	 */
	public void afterSubstParamChange(Tree tree, SubstitutionModel model, int ind, boolean accepted) {}
	
	/**
	 * Allows access to the Mcmc class. Generally not recommended.
	 */
	protected Mcmc getMcmc() {
		return manager.getMcmc();
	}
}
