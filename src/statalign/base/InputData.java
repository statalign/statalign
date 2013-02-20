package statalign.base;

import java.io.File;

import statalign.io.RawSequences;
import statalign.model.subst.SubstitutionModel;

public class InputData {

	/**
	 * The loaded sequences are stored in the RawSequences class
	 */
	public RawSequences seqs = new RawSequences();

	/**
	 * This class stores the parameters (number of burn-in cycles, number of
	 * steps, number of samplings) of the MCMC
	 */
	public MCMCPars pars = new MCMCPars(10000, 50000, 100, 1, 1, 100,
			new AutomateParamSettings());
	
	/**
	 * Are we using the parallel version?
	 */
	public boolean isParallel;

	/**
	 * The current substitution model that is used to analyse the sequences.
	 */
	public SubstitutionModel model;

	/**
	 * The path where output files will be created. By default this is the path of 
	 * the first sequence input file but can be configured by the user.
	 * @see #title
	 */
	public String outputPath;

	/**
	 * The file name (without path) that is chosen to be the base name for all output files,
	 * plugins will define their own output file names by adding extensions to this base name.
	 * By default, this is the name of the first sequence input file but can be configured
	 * by the user.
	 * @see #baseFilePath
	 */
	// TODO rename to baseFilename
	public String title;
	
	/**
	 * This integer stores the index of the current alignment type, as it is in
	 * <tt>alignmentTypes</tt>
	 */
	public int currentAlignmentType = 0;

	/**
	 * Fills both {@link #outputPath} and {@link #title}.
	 * @param file the name and path of the base file 
	 */
	public void setBaseFile(File file) {
		outputPath = file.getParent();
		title = file.getName();
	}
}
