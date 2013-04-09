package statalign.base.hmm;
import java.util.Arrays;

import statalign.base.Utils;

/**
 * An implementation of the abstract pair-HMM that emits characters into non-observable sequences
 * Used when a new alignment is proposed that aligns together two substrings via an ancestral sequence.
 * 
 * @author novak
 *
 */
public class HmmNonParam extends HmmSilent {

	//final double P = 0.9999999;
	//final double P = 0.99999; // original
	final double P = 0.99;
	final double Q = 0.6; // original
	//final double Q = 0.2;
	final int SILENT = 7;

	/* transition matrix for 3-sequence alignment HMM, st. 7 is silent, st. 0 is start, st. 6 is end*/
	private double transMatrix[][] = 
	{{0,P*P*P*P*P,P*P*P*P*(1-P),P*P*P*P*(1-P),1-P,P*(1-P),P-P*P*P-P*(1-P),P*P*P*(1-P)*(1-P)},
			{0,1-P+P*P*P*P*P*P,P*P*P*P*P*(1-P),P*P*P*P*P*(1-P),P*(1-P),P*P*(1-P),P*P*P*(1-P),P*P*P*P*(1-P)*(1-P)},
			{0,P*P*P*P*P*Q,1-Q+P*P*P*P*(1-P)*Q,P*P*P*P*(1-P)*Q,(1-P)*Q,P*(1-P)*Q,P*P*(1-P)*Q,P*P*P*(1-P)*(1-P)*Q},
			{0,P*P*P*P*P*Q,P*P*P*P*(1-P)*Q,1-Q+P*P*P*P*(1-P)*Q,(1-P)*Q,P*(1-P)*Q,P*P*(1-P)*Q,P*P*P*(1-P)*(1-P)*Q},
			{0,P*P*P*P*Q,P*P*P*(1-P)*Q,P*P*P*(1-P)*Q,1-Q,(1-P)*Q,P*(1-P)*Q,P*P*(1-P)*(1-P)*Q},
			{0,P*P*P*Q,P*P*(1-P)*Q,P*P*(1-P)*Q,0,1-Q,(1-P)*Q,P*(1-P)*(1-P)*Q},
			{0,0,0,0,0,0,0,0},
			{0,P*P*P*P*P*Q,P*P*P*P*(1-P)*Q,P*P*P*P*(1-P)*Q,(1-P)*Q,P*(1-P)*Q,P*P*(1-P)*Q,1-Q+P*P*P*(1-P)*(1-P)*Q}};

	/* reduced transition matrix, silent st. eliminated, st. 0 is start, st. 6 is end */
	private double redTransMatrix[][] = new double[7][7];

	/* states' emission descriptor: first dim. is parent/left child/right child, 2nd dim. is state, value is 0/1 */
	private int stateEmit[][] = {{0,1,1,1,0,0,0,1},
			{0,1,1,0,1,0,0,0},
			{0,1,0,1,0,1,0,0}};

	/* converts states' emission pattern as a binary number (p=4,l=2,r=1) to state # : e.g. 7->1, 6->2
       8->7 is virtual, so that end state has a pattern, too */
	private int emitPatt2State[] = {0,5,4,-1,7,3,2,1,6};

	/**
	 * Constructs a HMMSilent for Tree. Sets up the transition matrices.
	 */
	public HmmNonParam() {
		double sil2sil = 1-transMatrix[SILENT][SILENT];
		int i, j;

		for(i = 0; i <= 5; i++)
			for(j = 1; j < 7; j++)
				redTransMatrix[i][j] = Math.log(transMatrix[i][j]+transMatrix[i][SILENT]*transMatrix[SILENT][j]/sil2sil);
		for(i = 0; i < 7; i++) {
			redTransMatrix[i][0] = Utils.log0;
			redTransMatrix[6][i] = Utils.log0;
		}
		for(i = 0; i < 8; i++)
			for(j = 0; j < 8; j++)
				transMatrix[i][j] = Math.log(transMatrix[i][j]);
	}
	
	/**
	 * Returns the transition matrix.
	 * See the abstract class HmmSilent for more details.
	 */
	@Override
	public double[][] preCalcTransMatrix(double[][] transMatrix, double t1, double t2) {
		return this.transMatrix;
	}

	/**
	 * Returns the transition matrix.
	 * See the abstract class HmmSilent for more details.
	 */
	@Override
	public double[][] preCalcRedTransMatrix(double[][] redTransMatrix, double[][] transMatrix) {
		double[][] redTransMatrixCopy = new double[this.redTransMatrix.length][];
	
		for (int i = 0; i < this.redTransMatrix.length; i++) {
			redTransMatrixCopy[i] = this.redTransMatrix[i].clone();
		}
		return redTransMatrixCopy;
	}

	/**
	 * Returns the index of the silent state.
	 */
	@Override
	public int getSilent() {
		return SILENT;
	}
	
	/**
	 * Returns the index of the start state.
	 */
	@Override
	public int getStart() {
		return 0;
	}

	/**
	 * Returns the index of the end state.
	 */
	@Override
	public int getEnd() {
		return 6;
	}

	/**
	 * Returns an array specifying the emission pattern of each state. See class Hmm. 
	 */
	@Override
	public int[][] getStateEmit() {
		return stateEmit;
	}
	
	/**
	 * Returns a conversion array from emission patterns (coded as integers) into state
	 * indices. See class Hmm for details.
	 */
	@Override
	public int[] getEmitPatt2State() {
		return emitPatt2State;
	}

	/**
	 * For testing/debugging purposes.
	 * 
	 * @param args No argument is used.
	 */
	public static void main(String[] args) {
		HmmNonParam hmm = new HmmNonParam();
		double tm[][] = hmm.preCalcTransMatrix(null, 1, 1);
		double redtm[][] = hmm.preCalcRedTransMatrix(null, tm);

		double heat = 0.9;
		
		for (int i = 0; i < redtm.length; i++) {
			double tempSum = Utils.log0;
			for (int j = 0; j < redtm[i].length; j++) {
				redtm[i][j] = redtm[i][j] * heat;
				tempSum = Utils.logAdd(tempSum, redtm[i][j]);

			}
			if(tempSum != Double.NEGATIVE_INFINITY){
				for (int j = 0; j < redtm[i].length; j++) {
					redtm[i][j] = redtm[i][j] - tempSum;
				}
			}
		}

		boolean failed = false;
		for(int i = 0; i < 7; i++) {
			double sum = Utils.log0;
			for(int j = 0; j < 7; j++)
				sum = Utils.logAdd(sum, redtm[i][j]);
			sum = Math.exp(sum);
			System.out.println(Arrays.toString(redtm[i]));
			System.out.println("Row "+i+" sum: "+sum);
			if(Math.abs(sum-(i==6?0.0:1.0)) > 1e-5) {
				failed = true;
				for(int j = 0; j < 7; j++) 
					System.out.println(" "+i+"->"+j+" likelihood: "+Math.exp(redtm[i][j]));
			}
		}
		System.out.println(failed?"Test failed.":"Test passed.");
	}

}
