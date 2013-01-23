package statalign.model.ext;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import statalign.base.InputData;
import statalign.base.MainManager;
import statalign.base.Mcmc;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.hmm.Hmm;
import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.PluginParameters;

/**
 * Interface between StatAlign core and the {@link ModelExtension} plugins.
 * Not visible to plugins.
 * 
 * @author novak
 */
public class ModelExtManager {
	
	private MainManager mainMan;
	private Mcmc mcmc;
	
	private List<ModelExtension> pluginList;
	private List<ModelExtension> activeList;

	private int[] propWeights;
	/** Plugin selected for proposing parameter change */
	ModelExtension selectedPlugin;
	

	public ModelExtManager(MainManager mainMan) {
		this.mainMan = mainMan;
	}
	
	/**
	 * Discovers model extension plugins and initialises them.
	 */
	public void init(PluginParameters params) {
		pluginList = Utils.findPlugins(ModelExtension.class);
		for(ModelExtension plugin : pluginList) {
			plugin.setManager(this);
			plugin.init(params);
		}
	}
	
	/**
	 * Determines the list of active {@link ModelExtension} plugins.
	 * Calls {@link ModelExtension#initRun(statalign.base.InputData)} on each of them.
	 * @param inputData input data
	 */
	public void initRun(InputData inputData) throws IllegalArgumentException {
		activeList = new ArrayList<ModelExtension>();
		for(ModelExtension plugin : pluginList) {
			if(plugin.isActive()) {
				activeList.add(plugin);
				plugin.initRun(inputData);
			}
		}
		propWeights = new int[activeList.size()];
	}
	
	public void setMcmc(Mcmc mcmc) {
		this.mcmc = mcmc;
	}
	
	/**
	 * Allows access to the Mcmc class. Generally not recommended.
	 */
	public Mcmc getMcmc() {
		return mcmc;
	}
	
	/**
	 * Calls {@link ModelExtension#beforeSampling(Tree)} on each active plugin.
	 * @param tree current tree
	 */
	public void beforeSampling(Tree tree) {
		for(ModelExtension plugin : activeList)
			plugin.beforeSampling(tree);
	}
	
	public void afterSampling() {
		for(ModelExtension plugin : activeList)
			plugin.afterSampling();
	}
	
	/**
	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
	 * contributions from all model extension plugins to the log-likelihood of the tree.
	 * @param tree the current tree
	 * @return the total log-likelihood
	 */
	public double totalLogLike(Tree tree) {
		double ll = tree.getLogLike();
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeFactor(tree);
		return ll;
	}
	
	/**
	 * Calculates the log of total prior probability by including contributions from model extension plugins.
	 * @param tree the current tree
	 * @return the log of total prior probability 
	 */
	public double totalLogPrior(Tree tree) {
		double ll = tree.getLogPrior();
		for(ModelExtension plugin : activeList)
			ll += plugin.logPrior();
		return ll;
	}
	
	public int getParamChangeWeight() {
		int i = 0, total = 0;
		for(ModelExtension plugin : activeList)
			total += (propWeights[i++] = plugin.getParamChangeWeight());
		return total;
	}

	public void beforeModExtParamChange(Tree tree) {
		selectedPlugin = activeList.get(Utils.weightedChoose(propWeights));
		selectedPlugin.beforeModExtParamChange(tree, selectedPlugin);
		for(ModelExtension plugin : activeList)
			if(plugin != selectedPlugin)
				plugin.beforeModExtParamChange(tree, selectedPlugin);
	}
	
	public void proposeParamChange(Tree tree) {
		selectedPlugin.proposeParamChange(tree);
	}
	
	public boolean modExtParamChangeCallback(double logProposalRatio) {
		return mcmc.modExtParamChangeCallback(logProposalRatio);
	}
	
	public void modifyProposalWidths() {
		for(ModelExtension modExt : activeList)
			modExtParamModifyProposalWidths(modExt.proposalCounts, modExt.acceptanceCounts,
					modExt.proposalWidthControlVariables);
	}
	
	public void modExtParamModifyProposalWidths(
			int[] proposalCounts, int[] acceptanceCounts,
			double[] proposalWidthControlVariables) {
		
		for (int i=0; i<proposalCounts.length; i++) {
			if (proposalCounts[i] > Utils.MIN_SAMPLES_FOR_ACC_ESTIMATE) {
				double accRate = (proposalCounts[i] == 0) ? 0 : (double)acceptanceCounts[i]/(double)proposalCounts[i];
				//System.out.println("modExtParam["+i+"] acceptance rate = "+accRate);
				//System.out.println("var = "+proposalWidthControlVariables[i]);
				if (accRate < Utils.MIN_ACCEPTANCE) {
					proposalWidthControlVariables[i] *= Utils.SPAN_MULTIPLIER;
					proposalCounts[i] = 0;
					acceptanceCounts[i] = 0;
				}
				else if (accRate > Utils.MAX_ACCEPTANCE) {
					proposalWidthControlVariables[i] /= Utils.SPAN_MULTIPLIER;
					proposalCounts[i] = 0;
					acceptanceCounts[i] = 0;
				}
			}
		}
	}
	
//	/**
//	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
//	 * contributions from all model extension plugins to the log-likelihood of the tree.
//	 * @param tree the current tree
//	 * @return the total log-likelihood
//	 */
	public double logLikeModExtParamChange(Tree tree) {
		double ll = tree.getLogLike();
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeModExtParamChange(tree, selectedPlugin);
		return ll;
	}
	
	public void afterModExtParamChange(Tree tree, boolean accepted) {
		selectedPlugin.afterModExtParamChange(tree, selectedPlugin, accepted);
		for(ModelExtension plugin : activeList)
			if(plugin != selectedPlugin)
				plugin.afterModExtParamChange(tree, selectedPlugin, accepted);
	}
	
	public void beforeAlignChange(Tree tree, Vertex selectRoot) {
		for(ModelExtension plugin : activeList)
			plugin.beforeAlignChange(tree, selectRoot);
	}
	
//	/**
//	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
//	 * contributions from all model extension plugins to the log-likelihood of the tree.
//	 * @param tree the current tree
//	 * @return the total log-likelihood
//	 */
	public double logLikeAlignChange(Tree tree, Vertex selectRoot) {
		double ll = tree.getLogLike();
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeAlignChange(tree, selectRoot);
		return ll;
	}

	public void afterAlignChange(Tree tree, Vertex selectRoot, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterAlignChange(tree, selectRoot, accepted);
	}
	
	public void beforeTreeChange(Tree tree, Vertex nephew) {
		for(ModelExtension plugin : activeList)
			plugin.beforeTreeChange(tree, nephew);
	}
	
//	/**
//	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
//	 * contributions from all model extension plugins to the log-likelihood of the tree.
//	 * @param tree the current tree
//	 * @return the total log-likelihood
//	 */
	public double logLikeTreeChange(Tree tree, Vertex nephew) {
		double ll = tree.getLogLike();
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeTreeChange(tree, nephew);
		return ll;
	}
	
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterTreeChange(tree, nephew, accepted);
	}
	
	public void beforeEdgeLenChange(Tree tree, Vertex vertex) {
		for(ModelExtension plugin : activeList)
			plugin.beforeEdgeLenChange(tree, vertex);
	}
	
//	/**
//	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
//	 * contributions from all model extension plugins to the log-likelihood of the tree.
//	 * @param tree the current tree
//	 * @return the total log-likelihood
//	 */
	public double logLikeEdgeLenChange(Tree tree, Vertex vertex) {
		double ll = tree.getLogLike();
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeEdgeLenChange(tree, vertex);
		return ll;
	}
	
	public void afterEdgeLenChange(Tree tree, Vertex vertex, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterEdgeLenChange(tree, vertex, accepted);
	}
	
	public void beforeIndelParamChange(Tree tree, Hmm hmm, int ind) {
		for(ModelExtension plugin : activeList)
			plugin.beforeIndelParamChange(tree, hmm, ind);
	}
	
//	/**
//	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
//	 * contributions from all model extension plugins to the log-likelihood of the tree.
//	 * @param tree the current tree
//	 * @return the total log-likelihood
//	 */
	public double logLikeIndelParamChange(Tree tree, Hmm hmm, int ind) {
		double ll = tree.getLogLike();
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeIndelParamChange(tree, hmm, ind);
		return ll;
	}

	
	public void afterIndelParamChange(Tree tree, Hmm hmm, int ind, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterIndelParamChange(tree, hmm, ind, accepted);
	}
	
	public void beforeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {
		for(ModelExtension plugin : activeList)
			plugin.beforeSubstParamChange(tree, model, ind);
	}
	
//	/**
//	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
//	 * contributions from all model extension plugins to the log-likelihood of the tree.
//	 * @param tree the current tree
//	 * @return the total log-likelihood
//	 */
	public double logLikeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {
		double ll = tree.getLogLike();
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeSubstParamChange(tree, model, ind);
		return ll;
	}

	public void afterSubstParamChange(Tree tree, SubstitutionModel model, int ind, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterSubstParamChange(tree, model, ind, accepted);
	}

	
	/**
	 * Returns the (unmodifiable) list of recognised {@link ModelExtension} plugins.
	 */
	public List<ModelExtension> getPluginList() {
		return Collections.unmodifiableList(pluginList);
	}
	
}
