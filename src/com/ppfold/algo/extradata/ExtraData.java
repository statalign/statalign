package com.ppfold.algo.extradata;

import java.io.BufferedInputStream;
import java.util.List;

/**
 * Object class to contain extra data and processing functions.
 * The model for data interpretation has to be coded here. 
 * 
 * @author Z.Sukosd
 */


public interface ExtraData {
	
//	/**Returns the probability of observing the particular data 
//	 * given that the nucleotide in the input position is embedded in a stack. (Produced by
//	 * F->dFd). The position is given in ALIGNMENT coordinates so implementations of ExtraData
//	 * must return values relevant for the alignment, not individual sequences.  
//	 * @param position
//	 */
//	float getProbabilityGivenInnerPair(int position);
//	
//	/**Returns the probability of observing the particular data 
//	 * given that the nucleotide in the input position is involved in a basepair
//	 * at the end of a helix. (Produced by L->dFd) The position is given in ALIGNMENT
//	 * coordinates so implementations of ExtraData must return values relevant for the 
//	 * alignment, not individual sequences.  
//	 * @param position
//	 */
//	float getProbabilityGivenOuterPair(int position);
	
	/**Returns the probability of observing the particular data 
	 * given that the nucleotides in the input are paired with each other.
	 * The position is given in ALIGNMENT coordinates so implementations of ExtraData
	 * must return values relevant for the alignment, not individual sequences.  
	 * @param position1
	 * @param position2 
	 */
	
	float getProbabilityGivenOuterPaired(int position1, int position2);

	/**Returns the probability of observing the particular data 
	 * given that the nucleotide in the input position are paired with each other. 
	 * The position is given in ALIGNMENT coordinates so implementations of ExtraData
	 * must return values relevant for the alignment, not individual sequences.  
	 * @param position1
	 * @param position2 
	 */
	float getProbabilityGivenInnerPaired(int position1, int position2);

	/**Returns the probability of observing the particular data 
	 * given that the nucleotide in the input position is unpaired. 
	 * The position is given in ALIGNMENT coordinates so implementations of ExtraData
	 * must return values relevant for the alignment, not individual sequences.  
	 * @param position
	 */
	float getProbabilityGivenUnpaired(int position);
	

	/**Has to be able to transform its own dataset to remove the specified columns from it. 
	 * (These columns are gapped in the alignment) 
	 * @param leftoutcolumns
	 * @throws Exception 
	 */	
	void removeColumns(List<Integer> leftoutcolumns) throws Exception;

//	/**Has to be able to import data from a particular file. 
//	 * @param filename
//	 * @param sequencelength 
//	 * @throws Exception 
//	 */	
//	void importData(String filename, int sequencelength) throws Exception;
	
	/**Has to be able to import data from a BufferedInputStream. 
	 * @param stream
	 * @param sequencelength 
	 * @throws Exception 
	 */	
	void readData_toStream(BufferedInputStream stream, int sequencelength) throws Exception;
	
	/**Has to be able to transform the data for alignments. 
	 * @param gappedseq
	 */	
	void transformToAlignment(String gappedseq);

	/**Returns whether there is missing data data in the given position. 
	 * @param position
	 */	
	boolean isEmpty(int position);
	
	int getType();
	
	
}
