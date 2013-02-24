package statalign.postprocess;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import statalign.base.InputData;
import statalign.base.MainManager;
import statalign.base.Mcmc;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.base.Utils;
import statalign.distance.Pair;

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
	 * 
	 */
	public static boolean rnaMode = false;
	
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
	public void newSample(State state, int no, int total) {
		if(rnaMode) {
			for(Postprocess plugin : plugins) {
				plugin.newSample(state, no, total);
			}
		}
		
		else {
			for(Postprocess plugin : plugins) {
				if(!plugin.rnaAssociated) {
					plugin.newSample(state, no, total);
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
			for(Postprocess plugin : plugins)
				if(plugin.rnaAssociated)
					plugin.afterLastSample();
		} else {
			for(Postprocess plugin : plugins)
				if(!plugin.rnaAssociated)
					plugin.afterLastSample();
		}

	}
	
	public void initRun(InputData inputData) throws IOException {
		logFile = new FileWriter(new File(inputData.outputPath, inputData.title + ".log"));

		System.out.println();
		for (Postprocess p : plugins) {
			if (p.postprocessWrite && p.getFileExtension() != null) {
				String name = new File(inputData.outputPath, inputData.title + "." + p.getFileExtension()).getPath();
				System.out.println("Output file for " + p.getTabName()
						+ ": " + name);
				p.outputFile = new FileWriter(name);
			}
		}
	}
	
	public void finalizeRun() {
		try {
			logFile.close();
		} catch (IOException e) {
		}
		for (Postprocess p : plugins) {
			if (p.postprocessWrite && p.getFileExtension() != null) {
				try {
					p.outputFile.close();
				} catch (IOException e) {
				}
			}
		}

	}
	
	public List<Pair<String, String>> getFilesCreated() {
		List<Pair<String, String>> retList = new ArrayList<Pair<String, String>>();
		retList.add(new Pair<String, String>(mainManager.inputData.title+".log", "Log file with the MCMC samples"));
		for(Postprocess plugin : plugins) {
			List<String> files = plugin.getCreatedFileNames();
			List<String> dsc = plugin.getCreatedFileDescriptions();
			if(files == null) {
				String ext = plugin.getFileExtension();
				if(!plugin.postprocessWrite || ext == null)
					continue;
				files = Arrays.asList(mainManager.inputData.title + "." + ext);
			}
			for(int i = 0; i < files.size(); i++) {
				retList.add(new Pair<String, String>(files.get(i), 
						dsc != null && i < dsc.size() && dsc.get(i) != null ? 
						dsc.get(i) : "Output file for the '"+plugin.getTabName()+"' plugin"));
			}
		}
		return retList;
	}
	
}
