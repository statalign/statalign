package com.ppfold.algo;

import java.util.List;

/**
 * Job listener interface.
 * 
 * @author M.Vaerum
 */
public interface JobListener {
	/**
	 * Used in the case of CYKJob
	 */
	void jobFinished(JobResults results);

	/**
	 * Used in the case of PhyloJob
	 */
	void jobFinished(double[][] result);

	/**
	 * Used in the case of DataDistJob
	 */
	void jobFinished(List<ResultBundle> result);

}
