package com.ppfold.algo;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Definition of a job in the inside, outside and expectation algorithms.
 * 
 * @author Z.Sukosd
 */

public class CYKJob implements Serializable {
	private static final long serialVersionUID = 1L;

	// Inside algo: We only need to transport "L" for below and "S" for
	// diagbelow
	// Outside algo: Need to transport "S" and "F" for above and diagabove, "L"
	// for diagbelow
	// (different position in inside algo) and "S" for below (again different
	// position in inside algo)
	// therefore we need ResMatrix and not JobResults class

	// Expectation calc: We fill insideresultsS with BP probabilities,
	// below, diagbelow with expectation values, vertF, diagF also with
	// expectation values

	private List<ResMatrix> below = new ArrayList<ResMatrix>();
	private List<ResMatrix> diagbelow = new ArrayList<ResMatrix>();
	private List<ResMatrix> aboveS = new ArrayList<ResMatrix>();
	private List<ResMatrix> diagaboveS = new ArrayList<ResMatrix>();
	private List<ResMatrix> aboveF = new ArrayList<ResMatrix>();
	private List<ResMatrix> diagaboveF = new ArrayList<ResMatrix>();
	ResMatrix insideresultsL;
	ResMatrix insideresultsS;
	
	ResMatrix basepairs; //"normal" basepairs in case there's no difference between inner/outer, otherwise inner
	ResMatrix basepairs2; 
	private boolean diffBasepairs; //1 if inner and outer basepairs should be distinguished. 
	
	
	PointRes[] vertF;
	PointRes[] diagF;
	PointRes specialF;
	char specialntfirst;
	char specialntsecond;

	// following data only needed by outside algo
	PointRes[] vertL;
	PointRes[] diagL;
	PointRes specialL;

	int sectorid;

	int type; // 0 = inside; 1 = outside; 2 = expectation

	private int n;
	private int mini;
	private int minj;

	private int seqlength;
	private int fullseqlength;

	private double[][] prob; // rule probabilities

	public long getMemoryRequirement() {
		long resmatrixmemory = (long) (8 + (2 + (n * 12)) * 12 + 2 * 4 * n * n + 12);
		long memoryreq = 8 + 2 * 4; // serialversion + resmatrix references
		memoryreq += (insideresultsL != null ? 1 : 0) * resmatrixmemory; // insideResultL
		memoryreq += (insideresultsS != null ? 1 : 0) * resmatrixmemory; // insideResultS
		memoryreq += (basepairs != null ? 1 : 0) * resmatrixmemory; // basepairs
		memoryreq += below != null ? ((resmatrixmemory + 16) * below.size())
				: 0;
		memoryreq += diagbelow != null ? ((resmatrixmemory + 16) * diagbelow
				.size()) : 0;
		memoryreq += aboveS != null ? ((resmatrixmemory + 16) * aboveS.size())
				: 0;
		memoryreq += diagaboveS != null ? ((resmatrixmemory + 16) * diagaboveS
				.size()) : 0;
		memoryreq += aboveF != null ? ((resmatrixmemory + 16) * aboveF.size())
				: 0;
		memoryreq += diagaboveF != null ? ((resmatrixmemory + 16) * diagaboveF
				.size()) : 0;
		memoryreq += 2 * (4 + 4 + 12); // vertF/diagF - not checking if NULL but
		// marginal difference
		memoryreq += 7 * 4; // integers
		memoryreq += 12 + 12 * 2 + 6 * 8; // rule probabilities
		memoryreq *= 1.25; // add 25% extra space for safety
		return memoryreq;
	}

	public long getDataTransportRequirement() {
		long resmatrixmemory = (long) (8 + (2 + (n * 12)) * 12 + 2 * 4 * n * n + 12); // each
		// ResMatrix
		// uses
		// this
		// much
		// memory
		long memoryreq = 8 + 2 * 4; // serialversion + resmatrix references
		memoryreq += (insideresultsL != null ? 1 : 0) * resmatrixmemory; // insideResultL
		memoryreq += (insideresultsS != null ? 1 : 0) * resmatrixmemory; // insideResultS
		memoryreq += (basepairs != null ? 1 : 0) * resmatrixmemory; // basepairs
		memoryreq += below != null ? ((resmatrixmemory + 16) * below.size())
				: 0;
		memoryreq += diagbelow != null ? ((resmatrixmemory + 16) * diagbelow
				.size()) : 0;
		memoryreq += aboveS != null ? ((resmatrixmemory + 16) * aboveS.size())
				: 0;
		memoryreq += diagaboveS != null ? ((resmatrixmemory + 16) * diagaboveS
				.size()) : 0;
		memoryreq += aboveF != null ? ((resmatrixmemory + 16) * aboveF.size())
				: 0;
		memoryreq += diagaboveF != null ? ((resmatrixmemory + 16) * diagaboveF
				.size()) : 0;
		memoryreq += 2 * (4 + 4 + 12); // vertF/diagF - not checking if NULL but
		// marginal difference
		memoryreq += 7 * 4; // integers
		memoryreq += 12 + 12 * 2 + 6 * 8; // rule probabilities
		return memoryreq;
	}

	public long getResultTransportRequirement() {
		long resmatrixmemory = (long) (8 + (2 + (n * 12)) * 12 + 2 * 4 * n * n + 12);
		long memoryreq = 8 + 2 * 4; // serialversion + resmatrix references
		memoryreq += 3 * resmatrixmemory + 12; // results
		return memoryreq;
	}

	public long getExecutionTimeEstimate() {
		// all depend on job height. the scaling is NOT working!
		if (type == 0)
			// inside job
			//return (n * n * n) / (below.size() == 0 ? 24334 : 12167) + n
			//		* n * n * below.size() / 4000;
			return ((long)n*(long)n*(long)n)*(below.size()==0?1L:(below.size()*6L))/20000L;
		else if (type == 1)
			// outside job
			//return (n * n * n) / 2540 + n * n * n
			//		* (aboveS.size() + diagaboveS.size()) / 3200;
			return ((long)n*(long)n*(long)n)*(aboveS.size()+diagaboveS.size()+1L)/3000L;
		else
			// expectation job
			//return (n * n * n) / (below.size() == 0 ? 50334 : 24167) + n
			//		* n * n * below.size() / 13005;
			return ((long)n*(long)n*(long)n)*(below.size()==0?1L:(below.size()*6L))/50000L;
		
	}

	@Override
	public String toString() {
		if (type == 0) {
			return "inside " + sectorid;
		} else if (type == 1) {
			return "outside " + sectorid;
		} else {
			return "expectation " + sectorid;
		}
	}

	public int getSeqlength() {
		return seqlength;
	}

	public void setSeqlength(int seqlength) {
		this.seqlength = seqlength;
	}

	public int getFullseqlength() {
		return fullseqlength;
	}

	public void setFullseqlength(int fullseqlength) {
		this.fullseqlength = fullseqlength;
	}

	public int isType() {
		return type;
	}

	public void setType(int type) {
		this.type = type;
	}

	// Constructor
	public CYKJob(int sector, int no, int jobtype) {
		this.sectorid = sector;
		this.n = no;
		vertF = new PointRes[n];
		diagF = new PointRes[n];
		vertL = new PointRes[n];
		diagL = new PointRes[n];
		for (int i = 0; i < n; i++) {
			vertF[i] = new PointRes(0, 0);
			diagF[i] = new PointRes(0, 0);
		}
		specialF = new PointRes(0, 0);
		if (jobtype == 1) {
			for (int i = 0; i < n; i++) {
				vertL[i] = new PointRes(0, 0);
				diagL[i] = new PointRes(0, 0);
			}
			specialL = new PointRes(0, 0);
		}
		this.type = jobtype;
	}

	public int getSectorid() {
		return sectorid;
	}

	public int getMini() {
		return mini;
	}

	public int getMinj() {
		return minj;
	}

	public double getProb(int n, int m) {
		return prob[n][m];
	}

	public double[][] getProb() {
		return prob;
	}

	/**
	 * @return the below
	 */
	public List<ResMatrix> getBelow() {
		return below;
	}

	/**
	 * @param below
	 *            the below to set
	 */
	public void setBelow(List<ResMatrix> below) {
		this.below = below;
	}

	/**
	 * @return the diagbelow
	 */
	public List<ResMatrix> getDiagbelow() {
		return diagbelow;
	}

	/**
	 * @param diagbelow
	 *            the diagbelow to set
	 */
	public void setDiagbelow(List<ResMatrix> diagbelow) {
		this.diagbelow = diagbelow;
	}

	/**
	 * @return the aboveS
	 */
	public List<ResMatrix> getAboveS() {
		return aboveS;
	}

	/**
	 * @param aboveS
	 *            the aboveS to set
	 */
	public void setAboveS(List<ResMatrix> aboveS) {
		this.aboveS = aboveS;
	}

	/**
	 * @return the diagaboveS
	 */
	public List<ResMatrix> getDiagaboveS() {
		return diagaboveS;
	}

	/**
	 * @param diagaboveS
	 *            the diagaboveS to set
	 */
	public void setDiagaboveS(List<ResMatrix> diagaboveS) {
		this.diagaboveS = diagaboveS;
	}

	/**
	 * @return the aboveF
	 */
	public List<ResMatrix> getAboveF() {
		return aboveF;
	}

	/**
	 * @param aboveF
	 *            the aboveF to set
	 */
	public void setAboveF(List<ResMatrix> aboveF) {
		this.aboveF = aboveF;
	}

	/**
	 * @return the diagaboveF
	 */
	public List<ResMatrix> getDiagaboveF() {
		return diagaboveF;
	}

	/**
	 * @param diagaboveF
	 *            the diagaboveF to set
	 */
	public void setDiagaboveF(List<ResMatrix> diagaboveF) {
		this.diagaboveF = diagaboveF;
	}

	public void setData(int n, int mini, int minj, int seqlength,
			double[][] prob) {
		this.n = n;
		this.mini = mini;
		this.minj = minj;
		this.prob = prob;
		// this.probmatrix = probmatrix;
		this.seqlength = seqlength;
	}

	public void addBelowE(Sector sec) {
		below.add(sec.getExpectationMatrix());
	}

	public void addDiagBelowE(Sector sec) {
		diagbelow.add(sec.getExpectationMatrix());
	}

	public void addDiagBelowS(Sector sec, boolean jobtype) {
		if (!jobtype) {
			// inside job
			diagbelow.add(sec.getInsideMatrixS());
		} else {
			// outside job
			diagbelow.add(sec.getOutsideMatrixS());
		}
	}

	public void addBelowL(Sector sec, boolean jobtype) {
		if (!jobtype) {
			// inside job
			below.add(sec.getInsideMatrixL());
		} else {
			// outside job
			below.add(sec.getOutsideMatrixL());
		}
	}

	public void addBelowS(Sector sec, boolean jobtype) {
		// not used for inside
		below.add(sec.getInsideMatrixS());
	}

	public void addDiagBelowL(Sector sec, boolean jobtype) {
		// only used by outside
		diagbelow.add(sec.getInsideMatrixL());
	}

	public void addAboveS(Sector sec, boolean jobtype) {
		// only used by outside
		aboveS.add(sec.getOutsideMatrixS());

	}

	public void addDiagAboveS(Sector sec, boolean jobtype) {
		// only used by outside
		diagaboveS.add(sec.getOutsideMatrixS());
	}

	public void addAboveF(Sector sec, boolean jobtype) {
		// only used by outside
		aboveF.add(sec.getOutsideMatrixF());
	}

	public void addDiagAboveF(Sector sec, boolean jobtype) {
		// only used by outside
		diagaboveF.add(sec.getOutsideMatrixF());
	}

	public void setInsideVertF(Sector sec, int n, PointRes tmp) {
		vertF[n].copyFrom(sec.getMatrixInsideFTpVal(n, sec.dim - 1, tmp));
	}

	public void setInsideDiagF(Sector sec, int n, PointRes tmp) {
		diagF[n].copyFrom(sec.getMatrixInsideFTpVal(sec.dim - 1, n, tmp));
	}

	public void setExpectationVertE(Sector sec, int n, PointRes tmp) {
		vertF[n].copyFrom(sec.getMatrixExpectationTpVal(n, sec.dim - 1, tmp));
	}

	public void setExpectationDiagE(Sector sec, int n, PointRes tmp) {
		diagF[n].copyFrom(sec.getMatrixExpectationTpVal(sec.dim - 1, n, tmp));
	}

	public void setOutsideVertF(Sector sec, int n, PointRes tmp) {
		vertF[n].copyFrom(sec.getMatrixOutsideFTpVal(n, 0, tmp));
	}

	public void setOutsideDiagF(Sector sec, int n, PointRes tmp) {
		diagF[n].copyFrom(sec.getMatrixOutsideFTpVal(0, n, tmp));
	}

	public void setOutsideVertL(Sector sec, int n, PointRes tmp) {
		vertL[n].copyFrom(sec.getMatrixOutsideLTpVal(n, 0, tmp));
	}

	public void setOutsideDiagL(Sector sec, int n, PointRes tmp) {
		diagL[n].copyFrom(sec.getMatrixOutsideLTpVal(0, n, tmp));
	}

	public void setSpecialF(PointRes val) {
		specialF.copyFrom(val);
	}

	public PointRes getSpecialF() {
		return specialF;
	}

	public void setSpecialL(PointRes val) {
		specialL.copyFrom(val);
	}

	public PointRes getSpecialL() {
		return specialF;
	}

	public PointRes getVertF(int n) {
		return vertF[n];
	}

	public PointRes getDiagF(int n) {
		return diagF[n];
	}

	public PointRes getVertL(int n) {
		return vertL[n];
	}

	public PointRes getDiagL(int n) {
		return diagL[n];
	}

	public int getN() {
		return this.n;
	}

	public PointRes getBPProb(int s, int t) {
		return this.basepairs.getProb(s,t);
	}

	public void setDiffBasepairs(boolean diffBasepairs) {
		this.diffBasepairs = diffBasepairs;
	}

	public PointRes getinnerBPProb(int s, int t) {
		return this.basepairs.getProb(s,t);
	}
	
	public PointRes getouterBPProb(int s, int t) {
		if(diffBasepairs){
			//should return OUTER basepair prob
			return this.basepairs2.getProb(s,t);
		}
		else{
			//should return the only basepair prob
			return this.basepairs.getProb(s, t);
		}
	}
	
}
