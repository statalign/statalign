package statalign.base;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.SwingUtilities;

import mpi.MPI;
import statalign.base.mcmc.StatAlignMcmc;
import statalign.base.mcmc.StatAlignParallelMcmc;
import statalign.base.thread.MainThread;
import statalign.base.thread.StoppedException;
import statalign.io.DataManager;
import statalign.io.RawSequences;
import statalign.io.input.DataReader;
import statalign.model.ext.ModelExtManager;
import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.Postprocess;
import statalign.postprocess.PostprocessManager;
import statalign.postprocess.gui.InputGUI;
import statalign.ui.MainFrame;

/**
 * 
 * This is the central class in the information flow amongst classes.
 * 
 * This class also manages the main thread that runs the MCMC.
 * 
 * @author miklos, novak
 * 
 */
public class MainManager {

	/**
	 * Reference to all input data required to run StatAlign
	 */
	public InputData inputData = new InputData();
	
	/**
	 * This is the postprocess-manager of our program, which handles the
	 * postprocesses applied on the MCMC run.
	 */
	public PostprocessManager postProcMan;
	
	/**
	 * Data manager that handles all input files.
	 */
	public DataManager dataMan;
	
	/**
	 * Manager for model extension plugins
	 */
	public ModelExtManager modelExtMan;

	/**
	 * Array of substitution model classes that can be selected for an analysis.
	 */
	public Class<? extends SubstitutionModel>[] substModels;

	/**
	 * The main window of the program.
	 */
	public MainFrame frame;

	/**
	 * The graphical interface of the Input panel
	 */
	public InputGUI inputgui;

	/**
	 * Main (background) calculation thread of the application
	 */
	public MainThread thread;

	/**
	 * The full path of the input file from which we read the sequences.
	 */
	public String fullPath;
	
	/**
	 * After adding on extensions from ModelExtension objects.
	 */
	public String fullPathWithExtensions;

	/**
	 * Alignment formats in which StatAlign can generate output
	 * 
	 * Implemented formats currently are <tt>StatAlign</tt> (our own format),
	 * <tt>Clustal</tt>, <tt>Fasta</tt>, <tt>Phylip</tt>, <tt>Nexus</tt>
	 */
	public static String[] alignmentTypes = new String[] { "Fasta", "StatAlign",
			"Clustal",  "Phylip", "Nexus" };

	/**
	 * A trivial constructor that only sets <tt>MainFrame</tt> and substModels.
	 * 
	 * The MainFrame creates its MainManager, and the MainManager knows who is
	 * its owner MainFrame, so it can access GUIs on the MainFrame
	 * 
	 * @param frame
	 *            The owner of the MainManager, the main window of the graphical
	 *            interface.
	 */
	public MainManager(MainFrame frame, boolean parallel) {
		this.frame = frame;
		postProcMan = new PostprocessManager(this,parallel);
		dataMan = new DataManager();
		dataMan.init();
		modelExtMan = new ModelExtManager(this);
	}
	
	/**
	 * Initialises submanagers after the command line has been processed. This is separated from
	 * the constructor as some plugins require input from the command line.
	 * @param params plugin parameters 
	 */
	public void init(ArrayList<String> args) {
		modelExtMan.init(args);
		postProcMan.init(modelExtMan);
	}

	public void start() {
		start(1,1);
	}
	/**
	 * This function starts a new MCMC run.
	 * 
	 * It asks files for writing outputs, and it launches a <tt>MainThread</tt>
	 * that handles the MCMC run.
	 */
	public void start(int noOfProcesses, int rank) {
		
		try {
			String filenameExtension = modelExtMan.getFilenameExtension();
			if (noOfProcesses > 1) {
				filenameExtension += "chain"+rank;
			}
			if (!filenameExtension.isEmpty()) {
				filenameExtension += ".";
			}
			postProcMan.logFile = new FileWriter(fullPath + "." + filenameExtension + "log");					
			postProcMan.setBaseFileName(fullPath + "." + filenameExtension);
			postProcMan.createPluginFiles();					

			thread = new MainThread(this,rank,noOfProcesses);
			thread.start();
			
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Called when the MCMC thread terminates, signals end of the process back
	 * to MainFrame.
	 * @param errorCode -1: error 0: completed 1: stopped after sampling 2: stopped before sampling
	 */
	public void finished(final int errorCode, final Exception ex) {
		if (frame != null) {
			SwingUtilities.invokeLater(new Runnable() {
				@Override
				public void run() {
					frame.finished(errorCode, ex);
				}
			});
			
		}
	}
	
	public void deactivateRNA() {
		frame.deactivateRNA();
	}

}
