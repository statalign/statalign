package statalign.base;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import statalign.io.DataType;
import statalign.io.RawSequences;
import statalign.model.subst.SubstitutionModel;

public class InputData {

	/**
	 * The input sequences, optionally aligned
	 */
	public RawSequences seqs = new RawSequences();
	
	/**
	 * The initial tree, if given, which is optionally fixed for the whole analysis
	 */
	public Tree tree = null;

	/**
	 * Non-sequence type auxiliary data that has been loaded, can be used by model extension plugins
	 */
	public List<DataType> auxData = new ArrayList<DataType>();

	/**
	 * MCMC parameters such as the number of burn-in cycles, number of
	 * steps and sampling rate
	 */
	public MCMCPars pars = new MCMCPars(20000, 50000, 100, 1, 1, 100, 0, 0.05, 0.5);
	
	/**
	 * Specifies how the alignment of the input sequences in {@link #seqs} is used.<br><br>
	 *  0: it's ignored and sequences are re-aligned from scratch<br>
	 *  1: used as the initial alignment and sampled during MCMC<br>
	 *  2: fixed throughout the MCMC analysis (currently unsupported)
	 */
	public int useAlign = 0;

	/**
	 * Specifies how the given {@link #tree} is used.<br><br>
	 *  0: unused, the tree is constructed from scratch using NJ<br>
	 *  1: it is used as the initial tree and sampled during MCMC<br>
	 *  2: the topology is fixed throughout the MCMC analysis, edges are sampled<br>
	 *  3: both topology and edge lengths are fixed during the MCMC
	 */
	public int useTree = 0;
	
	/**
	 * Are we using the parallel version?
	 */
	public boolean isParallel;
	
	/** If <code>true</true> then logging information is also printed during the burnin. */
	public boolean doReportDuringBurnin = false;

	/**
	 * The current substitution model that is used to analyse the sequences.
	 */
	public SubstitutionModel model;
	
	/**
	* The file name (without path) that is chosen to be the base name for all output files,
	* plugins will define their own output file names by adding extensions to this base name.
	* By default, this is the name of the first sequence input file but can be configured
	* by the user.
	* @see #baseFilePath
	*/
	public String title;
	/**
	* The path where output files will be created. By default this is the path of 
	* the first sequence input file but can be configured by the user.
	* @see #title
	*/
	public String outputPath;

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
