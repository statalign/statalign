package statalign.base;

import java.util.ArrayList;
import java.util.List;

import statalign.io.DataType;
import statalign.io.RawSequences;
import statalign.model.subst.SubstitutionModel;

public class InputData {

	/**
	 * The loaded sequences are stored in the RawSequences class
	 */
	public RawSequences seqs = new RawSequences();

	/**
	 * Non-sequence type auxiliary data that has been loaded, can be used by model extension plugins
	 */
	public List<DataType> auxData = new ArrayList<DataType>();

	/**
	 * This class stores the parameters (number of burn-in cycles, number of
	 * steps, number of samplings) of the MCMC
	 */
	public MCMCPars pars = new MCMCPars(10000, 50000, 100, 1, 1, 100); // TODO:
																		// this.
	/**
	 * Are we using the parallel version?
	 */
	public boolean isParallel;

	/**
	 * The current substitution model that is used to analyse the sequences.
	 */
	public SubstitutionModel model;

	/**
	 * Title of the dataset (default: input filename without path)
	 */
	public String title;

	/**
	 * This integer stores the index of the current alignment type, as it is in
	 * <tt>alignmentTypes</tt>
	 */
	public int currentAlignmentType = 0;

}
