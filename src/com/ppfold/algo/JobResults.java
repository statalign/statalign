package com.ppfold.algo;

import java.io.Serializable;

/**
 * Contains the results of a CYKJob (3 ResMatrix objects)
 * 
 * @author Z.Sukosd
 * @see CYKJob
 */

public class JobResults implements Serializable {
	private static final long serialVersionUID = 1L;

	public ResMatrix L; // the result matrices of the job; expressing from L
	public ResMatrix S; // the result matrices of the job; expressing from S
	public ResMatrix F; // the result matrices of the job; expressing from F

	public JobResults() {
	}

	public JobResults(int n) {
		L = new ResMatrix(n);
		S = new ResMatrix(n);
		F = new ResMatrix(n);
	}

}
