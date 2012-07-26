package com.ppfold.algo;

import java.io.Serializable;

import com.sun.xml.internal.bind.v2.runtime.RuntimeUtil.ToStringAdapter;

/**
 * Contains matrices for storing inside, outside, expectation variables, and
 * methods for operations with them.
 * 
 * @author Z.Sukosd
 * @see JobResults
 * @see CYKJob
 */
public class ResMatrix implements Serializable {

	private static final long serialVersionUID = 1L;

	private float[][] fraction;
	private int[][] exponent;

	// dimensions of this resmatrix
	private int n;

	// Constructor
	public ResMatrix(int n) {
		fraction = new float[n][n];
		exponent = new int[n][n];
		this.n = n;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				fraction[i][j] = 0;
				exponent[i][j] = 0;
			}
		}
	}

	public void print() {
		for (int s = 0; s < n; s++) {
			for (int t = 0; t < n; t++) {
				System.out.print("(" + s + ", " + t + ")=" + fraction[s][t]
						+ "E" + exponent[s][t]);
			}
			System.out.println();
		}
	}

	public String printToString() {
		String result = "";
		for (int s = 0; s < n; s++) {
			for (int t = 0; t < n; t++) {
				result = result.concat("(" + s + ", " + t + ")="
						+ fraction[s][t] + "E" + exponent[s][t]);
			}
			result.concat("\n");
		}
		return result;
	}

	public PointRes getProb(int n, int m) {
		return new PointRes(fraction[n][m], exponent[n][m]);
	}

	public PointRes fetchProb(int n, int m, PointRes tmp) {
		tmp.setFraction(fraction[n][m]);
		tmp.setExponent(exponent[n][m]);
		return tmp;
	}

	public void setProb(int n, int m, PointRes point) {
		this.fraction[n][m] = point.getFraction();
		this.exponent[n][m] = point.getExponent();
	}

	public void addToProb(int n, int m, PointRes point) {
		this.setProb(n, m, this.getProb(n, m).add(point));
	}

	public void addToProb(int n, int m, PointRes point, PointRes tmp) {
		this.setProb(n, m, this.fetchProb(n, m, tmp).add(point));
	}

	public void replaceWithMaxSum(ResMatrix E1, ResMatrix E2, PointRes tmp,
			PointRes tmp2) {
		for (int s = 0; s < this.n; s++) {
			for (int t = 0; t < this.n; t++) {
				for (int ss = 0; ss < this.n; ss++) {
					tmp2.copyFrom(E1.fetchProb(s, ss, tmp));
					tmp2.add(E2.fetchProb(this.n - ss - 1, t, tmp));
					if (this.fetchProb(s, t, tmp).isLessThan(tmp2)) {
						this.setProb(s, t, tmp2);
					}

				}
			}
		}
	}

	public void incrementLWithDirectProduct(ResMatrix A, ResMatrix B,
			double prob, PointRes tmp, PointRes tmp2) {
		for (int s = this.n - 1; s >= 0; s--) {
			for (int t = this.n - 1; t >= 0; t--) {
				for (int k = 0; k <= this.n - 1; k++) {
					tmp.copyFrom(A.fetchProb(s, k, tmp2));
					tmp.multiply(B.fetchProb(this.n - (t + 1), k, tmp2), prob);
					this.addToProb(s, t, tmp, tmp2);
				}
			}
		}
	}

	public void incrementSWithDirectProduct(ResMatrix A, ResMatrix B,
			double prob, PointRes tmp, PointRes tmp2) {
		for (int s = this.n - 1; s >= 0; s--) {
			for (int t = this.n - 1; t >= 0; t--) {
				for (int k = 0; k <= this.n - 1; k++) {
					tmp.copyFrom(A.fetchProb(k, t, tmp2));
					tmp.multiply(B.fetchProb(k, this.n - (s + 1), tmp2), prob);
					this.addToProb(s, t, tmp, tmp2);
				}
			}
		}
	}

	public void incrementWithLS(ResMatrix L, ResMatrix S, double prob,
			PointRes tmp, PointRes tmp2) {
		for (int s = 0; s < this.n; s++) {
			for (int t = 0; t < this.n; t++) {
				for (int ss = 0; ss < this.n; ss++) {
					tmp.copyFrom(L.fetchProb(s, ss, tmp2));
					tmp.multiply(S.fetchProb(this.n - ss - 1, t, tmp2), prob);
					this.addToProb(s, t, tmp, tmp2);
				}
			}
		}
	}

	public String toString(){
		String matrixString = "";
		for (int s = 0; s < n; s++) {
			for (int t = 0; t < n; t++) {
				matrixString += fraction[s][t]	+ "E" + exponent[s][t] + "\t";
			}
			matrixString += "\n";
			
		}
		return matrixString;
	}

}
