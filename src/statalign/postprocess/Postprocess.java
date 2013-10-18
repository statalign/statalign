package statalign.postprocess;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JPanel;

import statalign.base.InputData;
import statalign.base.Mcmc;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.mcmc.McmcModule;
import statalign.model.ext.ModelExtManager;
import statalign.model.ext.ModelExtension;

/**
 * Common ancestor for post-processing plugin classes.
 * 
 * Comments for plugin developers:
 *  <ul><li> Implement the abstract class. You can use the boolean variables
 *           to tell the main program what your plugin is supposed to do.
 *      <li> Via the mcmc variable, you can reach the tree and eventually all
 *           its public variables, so you can collect information about the current state
 *           of the Markov chain
 *      <li> If you would like to make a GUI for this postprocess plugin, return with its
 *           JPanel in the getJPanel() function
 *      <li> The implemented subclass must handle the information gathered and needed 
 *           for the postprocessing. Neither the Mcmc class, nor the PostprocessManager class
 *           collects information about the state of the Markov chain!
 *  </ul>              
 *    
 */
public abstract class Postprocess {

	/** TODO: REMOVE!!!!!!!!!! */
	public Mcmc mcmc;
	
	PostprocessManager postprocessManager;
	
    /**
     * True, if the plugin implements the {@link Postprocess#getToolBarItems()} function.
     * Defaults to false.
     */
    public boolean hasToolBar = false;

	/**
	 * True if plugin is selected in the menu (and thus a tab is created for the plugin in
	 * the main window that can be used to allow the user to change settings before MCMC
	 * start and to show runtime information afterwards)
	 */
	public boolean selected = true;
	
	/**
	 * True if it <u>can</u> generate a GUI. Not used in the current version, it is
	 * for further development if one wants to switch on and off the GUIs.
	 */
	public boolean screenable = false;

	/**
	 * True if a GUI must be shown to the user.
	 */
	public boolean show = true;

	/**
	 * True if plugin is active (must produce its output) either because other plugins depend
	 * on it or because it is selected
	 */
	public boolean active = false;
	
	/**
	 * True if this class <u>can</u> generate an output
	 */
	public boolean outputable = false;

	/**
	 * True if this class <u>can</u> generate a postprocess file
	 */
	public boolean postprocessable = false;
	
	/**
	 * True if it writes into the log file
	 */
	public boolean sampling;
	
	/**
	 * True if it writes a postprocess file
	 */
	public boolean postprocessWrite;
	
	/**
	 * True if it should show only when running an RNA/DNA file
	 */
	public boolean rnaAssociated;
	
	
	
	/**
	 * This string tells the alignment type in which alignment must be presented
	 */
	public String alignmentType;
	
	/**
	 * This is the logfile writer that is written during the running and gets information from 
	 * all postprocesses.
	 */
	public FileWriter file;
	
	/**
	 * This is the output file writer, that is written by a specific postprocess.
	 */
	public FileWriter outputFile;
	

	/**
	 * Static object holding parameters, which are visible to all plugins.
	 */
	public static PluginParameters pluginParameters = new PluginParameters();


	/**
	 * Called by {@link PostprocessManager} after setting the {@link #show} field to allow initialisation,
	 * possibly involving GUI (toolbar etc.) 
	 */
//	public void init() {
//	}

	public void init(ModelExtManager modelExtMan) {
	}
	

	/**
	 * 
	 * @return Returns with the name that will appear on the label of the tabulated panel.
	 */
	public abstract String getTabName();
		
	/**
	 * 
	 * @return Returns with the icon that will appear on the label of the tabulated panel.
	 */
	public abstract Icon getIcon();
		
	/**
	 * 
	 * @return Returns with the panel of the GUI
	 */
	public abstract JPanel getJPanel();
	//public abstract Component getComponent();
	
	/**
	 * Reinitializes panel to accommodate a new GUI
	 */
	public void reloadPanel() {
	}
	
	
	/**
	 * Returns with the tip information (shown when the mouse cursor is moved over the 
	 * label of the tabulated panel)
	 */
	public abstract String getTip();

    /**
     * Specifies the order of the tab for this plugin in the GUI.
     * <p>
     * E.g. we probably always want the "sequence input"-tab to appear first and
     * therefore it returns a 1.0d here. The "current alignment"-tab returns 2.0d here
     * and so if we want to make a tab get ordered between the "sequence input" and the
     * "current alignment"-tabs we just make that plugin return 1.5d or something similar
     * here. By default tabs get ordered last.
     * @return the tab order of the tab associated with this plugin, if running in a GUI.
     */
    public double getTabOrder() {
        return Double.MAX_VALUE;
    }

	/**
	 * Returns default file extension that is to appended to the input
	 * file name to get the file this plugin is writing into.
	 */
	public String getFileExtension() {
		return null;
	}
	
	/**
	 * Override this and return an array of full-qualified class names of the plugins
	 * this plugin depends on.
	 */
	public String[] getDependences() {
		return null;
	}
	
	/**
	 * Override this to get access to instances of the plugins your plugin depends on.
	 * 
	 * This function will be called by the PostprocessManager during its initialisation.
	 * 
	 * @param plugins  reference to Postprocess objects in the order they are specified
	 *                 in getDependences() or null if it returns null
	 */
	public void refToDependences(Postprocess[] plugins) {
	}
	
	/**
	 * Called before MCMC start. This is the first time you can use PostprocessManager.mcmc
	 * to access internal data structure
	 */
	public void beforeFirstSample(InputData inputData) {
	}
	
	/**
	 * Called after a new step is made. A typical run of MCMC takes hundred thousands of
	 * steps, override this function only if it takes a negligible amount of time and does not
	 * use too much memory. We use, for example, in drawing the loglikelihood trace.
	 */
	public void newStep(McmcStep mcmcStep) {
	}
	
	/**
	 * Allows peeking into the Markov chain before actual sampling begins. Frequency is
	 * determined by the MCMC sampling parameter set by the user, just like for
	 * {@link #newSample(State, int, int)}. Unlike that method, this one is called both
	 * during and after the burn-in period.
	 * 
	 * @param state A {@link State} object representing the current state of the chain
	 */
	public void newPeek(State state) {
	}
	
	/**
	 * 
	 * This function is called when we sample from the Markov chain. Frequency is determined
	 * by the MCMC sampling parameter set by the user. Unlike {@link #newPeek(State state)},
     * this method is only called after the burn-in period.
	 * 
	 * @param state A {@link State} object representing the current state of the chain
	 * @param no The number of the current sample
	 * @param total The number of the total samples
	 */
	public void newSample(State state, int no, int total) {
	}
	public void newSample(McmcModule coreModel, State state, int no, int total) {		
		newSample(state,no,total);
	}
	/**
	 * This function switches on or off the sampling mode.
	 * @param enabled Set it true if you need samples.
	 */
	public abstract void setSampling(boolean enabled);
		
	
	/**
	 * This function is called after the MCMC runs.
	 */
	public void afterLastSample() {
	}

    /**
     * Returns the toolbar that this plugin uses.
     */
    public ArrayList<JComponent> getToolBarItems() {
        return null;
    }
    
    public List<ModelExtension> getModExtPlugins() {
    	return mcmc.modelExtMan.getPluginList();
    }
	
    public void setSelected(boolean selected) {
    	this.selected = selected;
    	postprocessManager.updateSelectedList();
    }

}
