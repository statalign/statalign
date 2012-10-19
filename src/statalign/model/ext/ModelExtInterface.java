package statalign.model.ext;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import statalign.base.InputData;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.hmm.Hmm;
import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.PluginParameters;

/**
 * Interface between StatAlign core and the {@link ModelExtension} plugins.
 * Not visible to plugins, utility methods are in {@link ModelExtManager}.
 * 
 * @author novak
 */
public class ModelExtInterface {
	
//	private MainManager mainMan;
	public ModelExtManager modelMan;
	
	private List<ModelExtension> pluginList;
	private List<ModelExtension> activeList;

	private int[] propWeights;
	
	public ModelExtInterface() {//MainManager mainMan) {
//		this.mainMan = mainMan;
		modelMan = new ModelExtManager(this);
	}
	
	/**
	 * Discovers model extension plugins and initialises them.
	 */
	public void init(PluginParameters params) {
		pluginList = Utils.findPlugins(ModelExtension.class);
		for(ModelExtension plugin : pluginList)
			plugin.init(modelMan, params);
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
	public double getTotalLogLike(Tree tree) {
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
	public double getTotalLogPrior(Tree tree) {
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
	
	public double proposeParamChange(Tree tree) {
		int ind = Utils.weightedChoose(propWeights);
		return activeList.get(ind).proposeParamChange(tree);
	}

	public void beforeAlignChange(Tree tree, Vertex selectRoot) {
		for(ModelExtension plugin : activeList)
			plugin.beforeAlignChange(tree, selectRoot);
	}
	
	public void afterAlignChange(Tree tree, Vertex selectRoot, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterAlignChange(tree, selectRoot, accepted);
	}
	
	public void beforeTreeChange(Tree tree, Vertex nephew) {
		for(ModelExtension plugin : activeList)
			plugin.beforeTreeChange(tree, nephew);
	}
	
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterTreeChange(tree, nephew, accepted);
	}
	
	public void beforeEdgeLenChange(Tree tree, Vertex vertex) {
		for(ModelExtension plugin : activeList)
			plugin.beforeEdgeLenChange(tree, vertex);
	}
	
	public void afterEdgeLenChange(Tree tree, Vertex vertex, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterEdgeLenChange(tree, vertex, accepted);
	}
	
	public void beforeIndelParamChange(Tree tree, Hmm hmm, int ind) {
		for(ModelExtension plugin : activeList)
			plugin.beforeIndelParamChange(tree, hmm, ind);
	}
	
	public void afterIndelParamChange(Tree tree, Hmm hmm, int ind, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterIndelParamChange(tree, hmm, ind, accepted);
	}
	
	public void beforeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {
		for(ModelExtension plugin : activeList)
			plugin.beforeSubstParamChange(tree, model, ind);
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
