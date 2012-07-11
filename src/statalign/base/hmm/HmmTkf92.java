package statalign.base.hmm;

import statalign.base.Utils;

/**
 * This class implements a simplified version of the TKF92 pair-HMM.
 * 
 * @author novak
 *
 */
public class HmmTkf92 extends Hmm2 {

	private double _B, _T, _R, _L, _M;
	
	private double HMM2_11() { return _R+(1-_R)*(1-_L*_B)*(_L/_M)*Math.exp(-_M*_T); }
	private double HMM2_12() { return (1-_R)*(1-_L*_B)*(_L/_M)*(1-Math.exp(-_M*_T)); }
	private double HMM2_13() { return (1-_R)*_L*_B; }
	private double HMM2_21() { return (1-_R)*(_L*_B/(1-Math.exp(-_M*_T)))*Math.exp(-_M*_T); }
	private double HMM2_22() { return _R+(1-_R)*_L*_B; }
	private double HMM2_23() { return (1-_R)*(1-_M*_B/(1-Math.exp(-_M*_T))); }
	private double HMM2_31() { return (1-_R)*(1-_L*_B)*_L/_M*Math.exp(-_M*_T); }
	private double HMM2_32() { return (1-_R)*(1-_L*_B)*_L/_M*(1-Math.exp(-_M*_T)); }
	private double HMM2_33() { return _R+(1-_R)*_L*_B; }
	private double HMM2_S1() { return (1-_L*_B)*_L/_M*Math.exp(-_M*_T); }
	private double HMM2_S2() { return (1-_L*_B)*_L/_M*(1-Math.exp(-_M*_T)); }
	private double HMM2_S3() { return _L*_B; }
	private double HMM2_SE() { return (1-_L*_B)*(1-_L/_M); }
	private double HMM2_1E() { return (1-_R)*(1-_L*_B)*(1-_L/_M); }
	private double HMM2_2E() { return (1-_R)*(_M*_B/(1-Math.exp(-_M*_T)))*(1-_L/_M); }
	private double HMM2_3E() { return (1-_R)*(1-_L*_B)*(1-_L/_M); }
	
	/* Emission (pattern) of the states: columns are states, 1st row: parent, 2nd row: child */
	private int stateEmit[][] = {{0,1,1,0,0},
															 {0,1,0,1,0}};
	
	/* converts states' emission pattern as a binary number (p=2,ch=1) to state # : e.g. 3->1, 2->2
	    4->4 is virtual, so that end state has a pattern, too */
	private int emitPatt2State[] = {0,3,2,1,4};

	/**
	 * This constructor creates a TKF92 pair-HMM for the tree
	 * @param defParams The default parameters of the model. Currently it is r = 0.2, lambda = 0.009,
	 * mu = 0.011.
	 */
	public HmmTkf92(double defParams[]) {
		if(defParams != null) {
			params = new double[defParams.length];
			System.arraycopy(defParams, 0, params, 0, defParams.length);
		} else {
			params = new double[] { 0.5, 0.009, 0.011 };			// TKF92 parameters: r, lambda, mu
			// !!! param init?
		}
	}
	
	/**
	 * Returns an array specifying the emission pattern of each state. See class Hmm. 
	 */
	public int[][] getStateEmit() {
		return stateEmit;
	}
	
	/**
	 * Returns a conversion array from emission patterns (coded as integers) into state
	 * indices. See class Hmm for details.
	 */
	public int[] getEmitPatt2State() {
		return emitPatt2State;
	}
	
	/**
	 * Returns the index of the start state.
	 */
	public int getStart() {
		return 0;
	}
	
	/**
	 * Returns the index of the end state.
	 */
	public int getEnd() {
		return 4;
	}

	/**
	 * Calculates a transition matrix given an edge length.
	 * See the abstract class Hmm2 for details.
	 */
	public double[][] preCalcTransMatrix(double[][] transMatrix, double t) {
		if(transMatrix == null)
			transMatrix = new double[5][5];		// TKF92 has 5 states including start (st. 0) & end (st. 4)

		_R = params[0]; _L = params[1]; _M = params[2];
		_T = t; _B = Math.exp((_L-_M)*_T); _B = (1-_B)/(_M-_L*_B);
		transMatrix[1][1] = Math.log(HMM2_11());
		transMatrix[1][2] = Math.log(HMM2_12());
		transMatrix[1][3] = Math.log(HMM2_13());
		transMatrix[2][1] = Math.log(HMM2_21());
		transMatrix[2][2] = Math.log(HMM2_22());
		transMatrix[2][3] = Math.log(HMM2_23());
		transMatrix[3][1] = Math.log(HMM2_31());
		transMatrix[3][2] = Math.log(HMM2_32());
		transMatrix[3][3] = Math.log(HMM2_33());
		transMatrix[0][1] = Math.log(HMM2_S1());
		transMatrix[0][2] = Math.log(HMM2_S2());
		transMatrix[0][3] = Math.log(HMM2_S3());
		transMatrix[0][4] = Math.log(HMM2_SE());
		transMatrix[1][4] = Math.log(HMM2_1E());
		transMatrix[2][4] = Math.log(HMM2_2E());
		transMatrix[3][4] = Math.log(HMM2_3E());
		for(int i = 0; i < 5; i++) {
			transMatrix[i][0] = Utils.log0;
			transMatrix[4][i] = Utils.log0;
		}
		return transMatrix;
	}

	/**
	 * For testing purposes.
	 * @param args No argument is used.
	 */
	public static void main(String args[]) {
		HmmTkf92 hmm = new HmmTkf92(null);
		double tm[][] = hmm.preCalcTransMatrix(null, 1);
		boolean failed = false;
		for(int i = 0; i < 5; i++) {
			double sum = Utils.log0;
			for(int j = 0; j < 5; j++)
				sum = Utils.logAdd(sum, tm[i][j]);
			sum = Math.exp(sum);
			System.out.println("Row "+i+" sum: "+sum);
			if(Math.abs(sum-(i==4?0.0:1.0)) > 1e-5) {
				failed = true;
				for(int j = 0; j < 5; j++) 
					System.out.println(" "+i+"->"+j+" likelihood: "+Math.exp(tm[i][j]));
			}
		}
		System.out.println(failed?"Test failed.":"Test passed.");
	}

}
