package statalign.model.ext;

import java.io.File;
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
import statalign.io.DataType;
import statalign.mcmc.McmcMove;
import statalign.model.subst.SubstitutionModel;

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
	
	private String filenameExtensionBase = "";

	public void addToFilenameExtension(String s) {
		if (filenameExtensionBase.isEmpty()) {
			filenameExtensionBase += s;
		}
		else {
			filenameExtensionBase += "_"+s;
		}
	}
	public String getFilenameExtension() {
		return filenameExtensionBase;
	}
	
	private int[] propWeights;
	/** Plugin selected for proposing parameter change */
	ModelExtension selectedPlugin;
	

	public ModelExtManager(MainManager mainMan) {
		this.mainMan = mainMan;
	}
	
	/**
	 * Discovers model extension plugins and initialises them.
	 */
	public void init(ArrayList<String> args) {
		pluginList = Utils.findPlugins(ModelExtension.class);
		for(ModelExtension plugin : pluginList) {
			plugin.setManager(this);
		}
		if (args != null) {
			ARG: for(int i = 0 ; i < args.size() ; ) {
				String [] pluginPair = args.get(i).split("\\[",2);
				for(ModelExtension plugin : pluginList) {
					if (plugin.getPluginID().equals(pluginPair[0])) {
						plugin.setActive(true);						
						if (pluginPair.length == 2) { // Then we have some parameters
							if (pluginPair[1].endsWith("]")) {
								String paramString = pluginPair[1].substring(0,pluginPair[1].length()-1);
								String [] paramList = paramString.split(",", -1);
								for (String param : paramList) {
									String [] paramPair = param.split("=", 2);
									if (paramPair.length == 2) {
										plugin.setParam(paramPair[0],paramPair[1]);
									}
									else {
										plugin.setParam(paramPair[0],true);
									}
								}
							}
							else {
								throw new IllegalArgumentException(
										"Plugin parameters must be specifed in the form\n-plugin:pluginName[par1=x,par2=y]\n");
							}
						}
						plugin.init();
						++i;
						continue ARG;
					}
				}
				String s = "";
				for(ModelExtension plugin : pluginList) {
					s += plugin.getPluginID()+" ("+(plugin.active?"active":"inactive")+")\n";
				}
				System.out.println("\nUnrecognised plugin: "+pluginPair[0]+"\n"
				+"Available plugins:\n"
				+s);				
				throw new IllegalArgumentException("Unrecognised plugin: "+pluginPair[0]);
			}
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
		for(ModelExtension plugin : pluginList) {
			plugin.setMcmc(mcmc);
		}
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
		//double ll = 0.0;
		for(ModelExtension plugin : activeList) {
			ll += plugin.logLikeFactor(tree);
		}
		return ll;
	}
	
	/**
	 * Calculates the log of total prior probability by including contributions from model extension plugins.
	 * @param tree the current tree
	 * @return the log of total prior probability 
	 */
	public double totalLogPrior(Tree tree) {
		double ll = tree.getLogPrior();
		//double ll = 0.0;
		for(ModelExtension plugin : activeList)
			ll += plugin.logPrior(tree);
		return ll;
	}
	
	public int getParamChangeWeight() {
		int i = 0, total = 0;
		for(ModelExtension plugin : activeList)
			total += (propWeights[i++] = plugin.getParamChangeWeight());
		return total;
	}
	
	public String getMcmcInfo() {
		String info = "";
		for(ModelExtension plugin : activeList) {
			info += "\n"+plugin.getPluginID()+"\n\n";
			info += plugin.getMcmcInfo();
		}
		return info;
	}
	public void beforeModExtParamChange(Tree tree) {
		selectedPlugin = activeList.get(Utils.weightedChoose(propWeights));
		selectedPlugin.beforeModExtParamChange(tree, selectedPlugin);
		for(ModelExtension plugin : activeList)
			if(plugin != selectedPlugin)
				plugin.beforeModExtParamChange(tree, selectedPlugin);
	}
	
	public boolean proposeParamChange(Tree tree) {
		return selectedPlugin.proposeParamChange(tree);
	}
	public boolean isParamChangeAccepted(double logProposalRatio,McmcMove m) {
		mcmc.coreModel.setLogLike(logLikeModExtParamChange(mcmc.tree));
		return mcmc.isParamChangeAccepted(logProposalRatio,m);
	}
	
	public void modifyProposalWidths() {
		for(ModelExtension modExt : activeList) {
			modExt.modifyProposalWidths();
		}
	}

	public void zeroAllMoveCounts() {
		for(ModelExtension modExt : activeList) {
			modExt.zeroAllMoveCounts();
		}
	}
	
//	public void resetAll() {
//		for(ModelExtension modExt : activeList) {
//			ModelExtension m = modExt.reset();
//			if (m != null) modExt = m;
//			m.init();
//		}
//	}
	/**
	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
	 * contributions from all model extension plugins to the log-likelihood of the tree.
	 * @param tree the current tree
	 * @return the total log-likelihood
	 */
	public double logLikeModExtParamChange(Tree tree) {
		double ll = tree.getLogLike();
		//double ll = 0.0;
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
	
	/**
	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
	 * contributions from all model extension plugins to the log-likelihood of the tree.
	 * @param tree the current tree
	 * @return the total log-likelihood
	 */
	public double logLikeAlignChange(Tree tree, Vertex selectRoot) {
		double ll = tree.getLogLike();
		//double ll = 0.0;
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
	
	/**
	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
	 * contributions from all model extension plugins to the log-likelihood of the tree.
	 * @param tree the current tree
	 * @return the total log-likelihood
	 */
	public double logLikeTreeChange(Tree tree, Vertex nephew) {
		double ll = tree.getLogLike();
		//double ll = 0.0;
		for(ModelExtension plugin : activeList) {
			ll += plugin.logLikeTreeChange(tree, nephew);	
		}
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
	
	/**
	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
	 * contributions from all model extension plugins to the log-likelihood of the tree.
	 * @param tree the current tree
	 * @return the total log-likelihood
	 */
	public double logLikeEdgeLenChange(Tree tree, Vertex vertex) {
		double ll = tree.getLogLike();
		//double ll = 0.0;
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeEdgeLenChange(tree, vertex);
		return ll;
	}
	
	public void afterEdgeLenChange(Tree tree, Vertex vertex, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterEdgeLenChange(tree, vertex, accepted);
	}
	
	public void beforeIndelParamChange(Tree tree, Hmm hmm, McmcMove m) {
		for(ModelExtension plugin : activeList)
			plugin.beforeIndelParamChange(tree, hmm, m);
	}
	
	/**
	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
	 * contributions from all model extension plugins to the log-likelihood of the tree.
	 * @param tree the current tree
	 * @return the total log-likelihood
	 */
	public double logLikeIndelParamChange(Tree tree, Hmm hmm, McmcMove m) {
		double ll = tree.getLogLike();
		//double ll = 0.0;
		for(ModelExtension plugin : activeList)
			ll += plugin.logLikeIndelParamChange(tree, hmm, m);
		return ll;
	}

	
	public void afterIndelParamChange(Tree tree, Hmm hmm, McmcMove m, boolean accepted) {
		for(ModelExtension plugin : activeList)
			plugin.afterIndelParamChange(tree, hmm, m, accepted);
	}
	
	public void beforeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {
		for(ModelExtension plugin : activeList)
			plugin.beforeSubstParamChange(tree, model, ind);
	}
	
	/**
	 * Calculates the total log-likelihood of the state by adding the log-likelihood factor
	 * contributions from all model extension plugins to the log-likelihood of the tree.
	 * @param tree the current tree
	 * @return the total log-likelihood
	 */
	public double logLikeSubstParamChange(Tree tree, SubstitutionModel model, int ind) {
		double ll = tree.getLogLike();
		//double ll = 0.0;
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
	public void afterFirstHalfBurnin() {
		for(ModelExtension plugin : activeList)
			plugin.afterFirstHalfBurnin();		
	}
	public void afterBurnin() {
		for(ModelExtension plugin : activeList)
			plugin.afterBurnin();		
	}
	public void incrementWeights() {
		for(ModelExtension plugin : activeList)
			plugin.incrementWeights();		
	}
	
	/**
	 * 
	 * @param aligned Vector indicating which characters are aligned to the current
	 * column in the subtrees below.
	 * @return Logarithm of emission probability for subtrees
	 */
	public double calcLogEm(int[] aligned) {
		double logEm = 0;
		for(ModelExtension plugin : activeList) {
			logEm += plugin.calcLogEm(aligned);
		}
		return logEm;
	}
	
	/**
	 * Called in GUI mode only when new data is added to the analysis.
	 * This is transferred to the plugins to allow them to respond.
	 */
	public void dataAdded(File file, DataType data) {
		for(ModelExtension plugin : pluginList) {
			plugin.dataAdded(file, data);
		}
	}		
}
