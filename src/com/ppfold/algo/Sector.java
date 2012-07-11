package com.ppfold.algo;

/**
 * Defines a sector and stores results on the submitter machine. CYKJobs are
 * created by Master when there is a need for them.
 * 
 * @author Z.Sukosd
 * @see SectorGenerator
 * @see Master
 */

public class Sector {

	public int dim; // dimensions of matrices in sector
	public Sector below; // the job directly below this one; L direction (job
	// depends on it)
	public Sector above; // the job directly above this one; L direction
	public Sector diagbelow; // the job diagonally below this one; S direction
	// (job depends on it)
	public Sector diagabove; // the job diagonally above this one; S direction
	public Sector next; // the one with one higher id

	// public int transferspaceInside; //the number of jobresults that need to
	// be transferred to calculate this sector
	// public int transferspaceOutside; //the number of jobresults that need to
	// be transferred to calculate this sector

	public int pos[]; // position of segment, defined by position of bottom left
	// corner

	private static final long serialVersionUID = 1L;

	private boolean isBeingInsideProcessed = false;
	private boolean isBeingOutsideProcessed = false;
	private boolean isBeingExpectationProcessed = false;

	private JobResults resultsInside; // the resultsInside matrices of the job
	// eg. resultsInside.L.p
	private JobResults resultsOutside;
	private JobResults resultsExpectation;
	private JobResults basepairs;
	private JobResults basepairs2;
	public boolean diffbp; //should inner and outer bp be differentiated? 
	
	public int sectorid;

	public int seqlength;
	public int fullseqlength; // length of the whole sequence (needed for

	// outside algo)

	public Sector() {
	}

	public void initializeBasepairs() {
		basepairs = new JobResults();
		basepairs.S = new ResMatrix(dim);
	}
	
	public void initializeBasepairs2() {
		basepairs2 = new JobResults();
		basepairs2.S = new ResMatrix(dim);
	}
	

	public void clearAllInside() {
		resultsInside = null;
	}

	public void clearAllOutside() {
		resultsOutside = null;
	}

	public void clearAllExpectation() {
		resultsExpectation = null;
	}

	public void clearAllBp() {
		basepairs = null;
	}

	public void setIsBeingInsideProcessed(boolean value) {
		isBeingInsideProcessed = value;
	}

	public boolean isBeingInsideProcessed() {
		return isBeingInsideProcessed;
	}

	public void setIsBeingOutsideProcessed(boolean value) {
		isBeingOutsideProcessed = value;
	}

	public boolean isBeingOutsideProcessed() {
		return isBeingOutsideProcessed;
	}

	public void setIsBeingExpectationProcessed(boolean value) {
		isBeingExpectationProcessed = value;
	}

	public boolean isBeingExpectationProcessed() {
		return isBeingExpectationProcessed;
	}

	public JobResults getInsideResults() {
		return resultsInside;
	}

	public JobResults getExpectationResults() {
		return resultsExpectation;
	}

	public JobResults getOutsideResults() {
		return resultsOutside;
	}

	public ResMatrix getExpectationMatrix() {
		return resultsExpectation.S;
	}

	public ResMatrix getBasePairs() {
		return basepairs.S;
	}

	public ResMatrix getBasePairs2() {
		return basepairs2.S;
	}
	
	public ResMatrix getInsideMatrixL() {
		return resultsInside.L;
	}

	public ResMatrix getInsideMatrixS() {
		return resultsInside.S;
	}

	public ResMatrix getInsideMatrixF() {
		return resultsInside.F;
	}

	public ResMatrix getOutsideMatrixF() {
		return resultsOutside.F;
	}

	public ResMatrix getOutsideMatrixS() {
		return resultsOutside.S;
	}

	public ResMatrix getOutsideMatrixL() {
		return resultsOutside.L;
	}

	public PointRes getMatrixInsideFTpVal(int n, int m) {
		return resultsInside.F.getProb(n, m);
	}

	public PointRes getMatrixInsideFTpVal(int n, int m, PointRes tmp) {
		tmp.copyFrom(resultsInside.F.fetchProb(n, m, tmp));
		return tmp;
		// return resultsInside.F.getTp(n, m);
	}

	public PointRes getMatrixOutsideFTpVal(int n, int m, PointRes tmp) {
		tmp.copyFrom(resultsOutside.F.fetchProb(n, m, tmp));
		return tmp;
		// return resultsOutside.F.getTp(n, m);
	}

	public PointRes getMatrixOutsideLTpVal(int n, int m, PointRes tmp) {
		tmp.copyFrom(resultsOutside.L.fetchProb(n, m, tmp));
		return tmp;
		// return resultsOutside.L.getTp(n, m);
	}

	public PointRes getMatrixExpectationTpVal(int n, int m, PointRes tmp) {
		tmp.copyFrom(resultsExpectation.S.fetchProb(n, m, tmp));
		return tmp;
	}

	public PointRes getMatrixBpVal(int n, int m, PointRes tmp) {
		tmp.copyFrom(basepairs.S.fetchProb(n, m, tmp));
		return tmp;
	}
	public PointRes getMatrixBpVal2(int n, int m, PointRes tmp) {
		tmp.copyFrom(basepairs2.S.fetchProb(n, m, tmp));
		return tmp;
	}

	
	
	public PointRes getMatrixExpectationTpVal(int n, int m) {
		return resultsExpectation.S.getProb(n, m);
	}

	public PointRes getMatrixBpVal(int n, int m) {
		return basepairs.S.getProb(n, m);
	}

	public PointRes getMatrixBpVal2(int n, int m) {
		return basepairs2.S.getProb(n, m);
	}
	
	public void setMatrixExpectationTpVal(int n, int m, PointRes val) {
		resultsExpectation.S.setProb(n, m, val);
	}

	public void setBasePairs(int n, int m, PointRes val) {
		basepairs.S.setProb(n, m, val);
	}
	
	public void setBasePairs2(int n, int m, PointRes val) {
		basepairs2.S.setProb(n, m, val);
	}

	public void setDependency(Sector n, Sector m) {
		below = n;
		diagbelow = m;
	}

	public void setNextJobs(Sector n, Sector m) {
		above = n;
		diagabove = m;
	}

	public void setPos(int i, int j) {
		pos = new int[2];
		pos[0] = i;
		pos[1] = j;
	}

	public int getMini() {
		return pos[0];
	}

	public int getMinj() {
		return pos[1];
	}

	public void setInsideResults(JobResults results) {

		this.resultsInside = results;
	}

	public void setOutsideResults(JobResults results) {

		this.resultsOutside = results;
	}

	public void setExpectationResults(JobResults results) {

		this.resultsExpectation = results;
	}

	public boolean readyForInsideProcessing() {
		if ((this.getInsideResults() == null && !isBeingInsideProcessed())
				&& (diagbelow == null || diagbelow.getInsideResults() != null)
				&& (below == null || below.getInsideResults() != null)) {
			return true;
		} else {
			return false;
		}
	}

	public boolean readyForOutsideProcessing() {
		if ((this.getOutsideResults() == null && !isBeingOutsideProcessed())
				&& (diagabove == null || diagabove.getOutsideResults() != null)
				&& (above == null || above.getOutsideResults() != null)) {
			return true;
		} else {
			return false;
		}
	}

	public boolean readyForExpectationProcessing() {
		if ((this.getExpectationResults() == null && !isBeingExpectationProcessed())
				&& (diagbelow == null || diagbelow.getExpectationResults() != null)
				&& (below == null || below.getExpectationResults() != null)) {
			return true;
		} else {
			return false;
		}
	}

	public PointRes topinsideresult(int seqlength) {
		return this.getInsideMatrixS().getProb(dim + pos[0] - 1,
				seqlength - 1 - pos[0] - pos[1]);
	}

}
