package statalign.postprocess;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import statalign.base.MainManager;
import statalign.base.Mcmc;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.base.Utils;
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
	public Postprocess[] plugins;
	
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
	 * This constructor recognizes the plugins
	 * @param mainManager The MainManager that manages the MCMC run.
	 */
	public PostprocessManager(MainManager mainManager) {
		this.mainManager = mainManager;
		String[] pluginNames = Utils.classesInPackage(Postprocess.class.getPackage().getName()+".plugins");
		HashMap<String,Integer> nameMap = new HashMap<String,Integer>();
		plugins = new Postprocess[pluginNames.length];
		for(int i = 0; i < pluginNames.length; i++) {
			nameMap.put(pluginNames[i], i);
			try {
				Class<?> cl = Class.forName(pluginNames[i]);
				if(!Postprocess.class.isAssignableFrom(cl))
					continue;
				plugins[i] = (Postprocess)cl.newInstance();
				plugins[i].selected = plugins[i].active = true;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		ArrayList<Postprocess> workingPlugins = new ArrayList<Postprocess>();
		for(Postprocess plugin : plugins)
			dependProb(plugin, nameMap, workingPlugins);

		plugins = workingPlugins.toArray(new Postprocess[workingPlugins.size()]);

		boolean show = mainManager.frame != null;
		for(Postprocess plugin : plugins){
			plugin.selected = true;
			if(plugin.screenable) {
				plugin.show = show;
			}
			plugin.init();
		}
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
				if(ind == null || dependProb(plugins[ind], nameMap, workingPlugins))
					return true;
				refs[i] = plugins[ind];
			}
		}
		plugin.refToDependences(refs);
		plugin.active = false;
		workingPlugins.add(plugin);
		return false;
	}
	
	/**
	 * This function invoked before the first sample.
	 * It opens information chanels towards postrocessing plug-ins.
	 */
	public void beforeFirstSample() {
		for(Postprocess plugin : plugins){
			plugin.mcmc = mcmc;
			plugin.file = logFile;
			plugin.alignmentType = MainManager.alignmentTypes[mainManager.inputData.currentAlignmentType];
			plugin.beforeFirstSample(mainManager.inputData);
		}
	}
	
	/**
	 * Calls the plug-ins at a new sample.
	 * @param no The number of the current sample
	 * @param total The total number of samples.
	 */
	public void newSample(State state, int no, int total) {
		for(Postprocess plugin : plugins)
			plugin.newSample(state, no, total);
	}
	
	/**
	 * Calls the plug-ins after an MCMC step.
	 */
	public void newStep(McmcStep step) {
		for(Postprocess plugin : plugins){
			//if(plugin.selected){
				plugin.newStep(step);
			//}
		}
	}
	
	/**
	 * Calls the plug-ins after an MCMC step.
	 */
	public void newPeek() {
		State state = mcmc.getState();
		for(Postprocess plugin : plugins){
			//if(plugin.selected){
				plugin.newPeek(state);
			//}
		}
	}
	
	/**
	 * It is called after the last sample of MCMC.
	 * It calls plug-ins to finalize their postprocessing activities.
	 */
	public void afterLastSample() {
		for(Postprocess plugin : plugins){
			if(plugin.postprocessWrite){
				plugin.afterLastSample();
			}
		}

		try {
			logFile.close();
		} catch (IOException e) {
			new ErrorMessage(mainManager.frame, " "+e.getLocalizedMessage(), true);
		}
	}
	
}
