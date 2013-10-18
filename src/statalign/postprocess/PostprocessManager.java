package statalign.postprocess;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import statalign.base.MainManager;
import statalign.base.Mcmc;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.base.Utils;
import statalign.mcmc.McmcModule;
import statalign.model.ext.ModelExtManager;
import statalign.ui.ErrorMessage;

/**
 * This class manages the postprocesses.
 * 
 * @author miklos, novak
 *
 */
public class PostprocessManager {
	
	/**
	 * The recognized plugins are in this array
	 */
	private List<Postprocess> plugins;
	
	/**
	 * This is the Mcmc that is analyzed
	 */
	public Mcmc mcmc;
	
	/**
	 * This is the log-file writer 
	 */
	public FileWriter logFile = null;
	
	/**
	 * This is the Main manager that manages the MCMC run.
	 */
	public MainManager mainManager;
	
	/**
	 * 
	 */
	public static boolean rnaMode = false;
	
	/**
	 * Static object holding parameters, which are visible to all plugins.
	 */
	public static PluginParameters pluginParameters = new PluginParameters();

	
	/**
	 * This constructor recognizes the plugins
	 * @param mainManager The MainManager that manages the MCMC run.
	 */
	public PostprocessManager(MainManager mainManager) {
		this.mainManager = mainManager;
		List<String> pluginNames = Utils.classesInPackage(Postprocess.class.getPackage().getName()+".plugins");
		HashMap<String,Integer> nameMap = new HashMap<String,Integer>();
		HashMap<Postprocess,Boolean> selectedMap = new HashMap<Postprocess, Boolean>();
		plugins = new ArrayList<Postprocess>(pluginNames.size());
		for(int i = 0; i < pluginNames.size(); i++)
			plugins.add(null);
		for(int i = 0; i < pluginNames.size(); i++) {
			nameMap.put(pluginNames.get(i), i);
			try {
				Class<?> cl = Class.forName(pluginNames.get(i));
				if(!Postprocess.class.isAssignableFrom(cl))
					continue;
				Postprocess plugin = (Postprocess)cl.newInstance();
				selectedMap.put(plugin, plugin.selected);
				plugin.selected = plugin.active = true;
				plugins.set(i, plugin);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		ArrayList<Postprocess> workingPlugins = new ArrayList<Postprocess>();
		for(Postprocess plugin : plugins)
			dependProb(plugin, nameMap, workingPlugins);

		plugins = workingPlugins;

		boolean show = mainManager.frame != null;
		for(Postprocess plugin : plugins){
			plugin.selected = selectedMap.get(plugin);
			if(plugin.screenable) {
				plugin.show = show;
			}
//			plugin.init();
		}
	}
	
	public List<Postprocess> getPlugins() {
		return Collections.unmodifiableList(plugins);
	}
	
	/**
	 * Finds dependency problems. Adds working plugins in dependency order to workingPlugins.
	 */
	private boolean dependProb(Postprocess plugin, HashMap<String,Integer> nameMap, ArrayList<Postprocess> workingPlugins) {
		if(plugin == null)
			return true;
		if(!plugin.selected)
			return plugin.active;
		plugin.selected = false;
		String[] deps = plugin.getDependences();
		Postprocess[] refs = null;
		if(deps != null) {
			refs = new Postprocess[deps.length];
			for(int i = 0; i < deps.length; i++) {
				Integer ind = nameMap.get(deps[i]);
				if(ind == null || dependProb(plugins.get(ind), nameMap, workingPlugins))
					return true;
				refs[i] = plugins.get(ind);
			}
		}
		plugin.refToDependences(refs);
		plugin.active = false;
		workingPlugins.add(plugin);
		return false;
	}
	
	public void init(ModelExtManager modelExtMan) {
		for(Postprocess plugin : plugins){
			plugin.postprocessManager = this;
			plugin.init(modelExtMan);
		}
	}
	/**
	 * This function invoked before the first sample.
	 * It opens information chanels towards postrocessing plug-ins.
	 */
	public void beforeFirstSample() {
		if(rnaMode) {
			for(Postprocess plugin : plugins){
				plugin.mcmc = mcmc;
				plugin.file = logFile;
				plugin.alignmentType = MainManager.alignmentTypes[mainManager.inputData.currentAlignmentType];
				plugin.beforeFirstSample(mainManager.inputData);
			}
		}
		
		else {
			for(Postprocess plugin : plugins) {
				if(!plugin.rnaAssociated) {
					plugin.mcmc = mcmc;
					plugin.file = logFile;
					plugin.alignmentType = MainManager.alignmentTypes[mainManager.inputData.currentAlignmentType];
					plugin.beforeFirstSample(mainManager.inputData);
				}
			}
		}
	}
	
	/**
	 * Calls the plug-ins at a new sample.
	 * @param no The number of the current sample
	 * @param total The total number of samples.
	 */
	public void newSample(McmcModule coreModel, State state, int no, int total) {
		if(rnaMode) {
			for(Postprocess plugin : plugins) {
				plugin.newSample(coreModel,state, no, total);
			}
		}
		
		else {
			for(Postprocess plugin : plugins) {
				if(!plugin.rnaAssociated) {
					plugin.newSample(coreModel,state, no, total);
				}
			}
		}
	}
	
	/**
	 * Calls the plug-ins after an MCMC step.
	 */
	public void newStep(McmcStep step) {
		if(rnaMode) {
			for(Postprocess plugin : plugins){
				//if(plugin.selected){
					plugin.newStep(step);
				//}
			}
		}
		
		else {
			for(Postprocess plugin : plugins){
				if(!plugin.rnaAssociated){
					plugin.newStep(step);
				}
			}
		}
	}
	
	public void flushAll() {
		if(rnaMode) {
			for(Postprocess plugin : plugins){
				if(plugin.selected && plugin.postprocessWrite){
					try { plugin.outputFile.flush(); }
					catch (IOException ioex) {
						ioex.printStackTrace();
					} 
				}
			}
		}
		
		else {
			for(Postprocess plugin : plugins){
				if(plugin.selected && plugin.postprocessWrite && 
						!plugin.rnaAssociated){
					try { plugin.outputFile.flush(); }
					catch (IOException ioex) {
						ioex.printStackTrace();
					}
				}
			}
		}		
	}
	/**
	 * Calls the plug-ins after an MCMC step.
	 */
	public void newPeek() {
		State state = mcmc.getState();
		if(rnaMode) {
			for(Postprocess plugin : plugins){
				//if(plugin.selected){
					plugin.newPeek(state);
				//}
			}
		}
		
		else {
			for(Postprocess plugin : plugins){
				if(!plugin.rnaAssociated){
					plugin.newPeek(state);
				}
			}
		}
	}
	
	/**
	 * It is called after the last sample of MCMC.
	 * It calls plug-ins to finalize their postprocessing activities.
	 */
	public void afterLastSample() {
		if(rnaMode) {
			for(Postprocess plugin : plugins) {
				if(plugin.rnaAssociated) {
					plugin.afterLastSample();
				}

			}
		}
		else {
			for(Postprocess plugin : plugins) {
				if(!plugin.rnaAssociated) {
					plugin.afterLastSample();
				}
			}
		}			
		try {
			logFile.close();
		} catch (IOException e) {
			new ErrorMessage(mainManager.frame, " "+e.getLocalizedMessage(), true);
		}
	}

	public void updateSelectedList() {
		if(mainManager.frame != null)
			mainManager.frame.updateTabs();
	}
	
}
