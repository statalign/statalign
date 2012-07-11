package com.ppfold.algo;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.List;

/**
 * Reads and stores the input parameters (matrices) of the algorithm
 * 
 * @author Z.Sukosd
 */

public class Parameters implements Serializable {
	
	private static final long serialVersionUID = -6960312420312527995L;

	private boolean entropycalc = false; //Change to enable or disable calculation of entropy
	
	// rules of the grammar
	private double[][] prob;

	// reduced group
	// probabilities
	private double[] pr;
	// substitution rates diagonalized 
	private double[][] rV;   
	private double[][] rD;  
	private double[][] rV1; 
	
	// single group
	// unpaired probabilities
	private double[] ps;
	// unpaired substitution rates diagonalized 
	private double[][] sV;   
	private double[][] sD;  
	private double[][] sV1; 
	

	// double group
	// paired probabilities
	private double[][] pd;
	// paired bases substitution rates diagonalized 
	private double[][] dV;
	private double[][] dD;
	private double[][] dV1;

	
	
	public Parameters(double[][] prob, double[] pr, double[] ps, double[][] pd,
			double[][] rV, double[][] rD, double[][] rV1,
			double[][] sV, double[][] sD, double[][] sV1, double[][] dV,
			double[][] dD, double[][] dV1) {
		this.prob = prob;
		this.ps = ps;
		this.pd = pd;
		this.pr = pr;
		this.rV = rV;
		this.rD = rD;
		this.rV1 = rV1;
		this.sV = sV;
		this.sD = sD;
		this.sV1 = sV1;
		this.dV = dV;
		this.dD = dD;
		this.dV1 = dV1;

	}
	
	public static Parameters readParam(List<String> lines) throws IOException {
		double[][] prob = new double[3][2];

		
		// reduced matrices 
		double[] pr = new double[4];
		double[][] rV = new double[4][4];
		double[][] rD = new double[4][4];
		double[][] rV1 = new double[4][4];
		
		// single group
		double[] ps = new double[4];
		double[][] sV = new double[4][4];
		double[][] sD = new double[4][4];
		double[][] sV1 = new double[4][4];

		// double group
		double[][] pd = new double[4][4];
		double[][] dV = new double[16][16];
		double[][] dD = new double[16][16];
		double[][] dV1 = new double[16][16];

		int i = 0;

		int state = 0;
		int substate = 0;

		for (String line : lines) {
			if (line.length() > 0 && line.charAt(0) == '>') {
				state = state + 1;
				substate = 0;
				i = 0;
			} else {
				if (line.trim().length() == 0) {
					substate += 1;
					i = 0;
				} else {
					double[] d = line2double(line);
					switch (state) {
					// Grammar
					case 1:
						for (int k = 0; k < d.length; k++) {
							prob[k][0] = d[k];
							prob[k][1] = 1 - d[k];
						}
						break;
					// Reduced group
					case 2:
						switch (substate) {
						// 4x4
						case 0:
							pr = d;
							break;
						// 4x4
						case 1:
							rV[i] = d;
							break;
						// 4x4
						case 2:
							rD[i] = d;
							break;
						// 4x4
						case 3:
							rV1[i] = d;
							break;
						default:
							break;
						}
						break;
					case 3:
						switch (substate) {
						// 4x4
						case 0:
							ps = d;
							break;
						// 4x4
						case 1:
							sV[i] = d;
							break;
						// 4x4
						case 2:
							sD[i] = d;
							break;
						// 4x4
						case 3:
							sV1[i] = d;
							break;
						default:
							break;
						}
						break;
					// Double group
					case 4:
						switch (substate) {
						// 4x4
						case 0:
							pd[i] = d;
							break;
						// 16x16
						case 1:
							dV[i] = d;
							break;
						// 16x16
						case 2:
							dD[i] = d;
							break;
						// 16x16
						case 3:
							dV1[i] = d;
							break;
						default:
							break;
						}
						break;
					default:
						break;
					}
					i += 1;
				}
			}
		}

		return new Parameters(prob, pr, ps, pd, rV, rD, rV1, sV, sD, sV1, dV, dD, dV1);

	}

	public static Parameters readParam(BufferedReader input) throws IOException {
		double[][] prob = new double[3][2];

		
		// reduced matrices 
		double[] pr = new double[4];
		double[][] rV = new double[4][4];
		double[][] rD = new double[4][4];
		double[][] rV1 = new double[4][4];
		
		// single group
		double[] ps = new double[4];
		double[][] sV = new double[4][4];
		double[][] sD = new double[4][4];
		double[][] sV1 = new double[4][4];

		// double group
		double[][] pd = new double[4][4];
		double[][] dV = new double[16][16];
		double[][] dD = new double[16][16];
		double[][] dV1 = new double[16][16];

		int i = 0;

		int state = 0;
		int substate = 0;
		String line = "";

		while ((line = input.readLine()) != null) {
			if (line.length() > 0 && line.charAt(0) == '>') {
				state = state + 1;
				substate = 0;
				i = 0;
			} else {
				if (line.trim().length() == 0) {
					substate += 1;
					i = 0;
				} else {
					double[] d = line2double(line);
					switch (state) {
					// Grammar
					case 1:
						for (int k = 0; k < d.length; k++) {
							prob[k][0] = d[k];
							prob[k][1] = 1 - d[k];
						}
						break;
					// Reduced group
					case 2:
						switch (substate) {
						// 4x4
						case 0:
							pr = d;
							break;
						// 4x4
						case 1:
							rV[i] = d;
							break;
						// 4x4
						case 2:
							rD[i] = d;
							break;
						// 4x4
						case 3:
							rV1[i] = d;
							break;
						default:
							break;
						}
						break;
					case 3:
						switch (substate) {
						// 4x4
						case 0:
							ps = d;
							break;
						// 4x4
						case 1:
							sV[i] = d;
							break;
						// 4x4
						case 2:
							sD[i] = d;
							break;
						// 4x4
						case 3:
							sV1[i] = d;
							break;
						default:
							break;
						}
						break;
					// Double group
					case 4:
						switch (substate) {
						// 4x4
						case 0:
							pd[i] = d;
							break;
						// 16x16
						case 1:
							dV[i] = d;
							break;
						// 16x16
						case 2:
							dD[i] = d;
							break;
						// 16x16
						case 3:
							dV1[i] = d;
							break;
						default:
							break;
						}
						break;
					default:
						break;
					}
					i += 1;
				}
			}
		}

		return new Parameters(prob, pr, ps, pd, rV, rD, rV1, sV, sD, sV1, dV, dD, dV1);

	}

	private static double[] line2double(String line) {
		if ("".equals(line))
			return new double[0];
		// ltrim
		while (line.charAt(0) == ' ')
			line = line.substring(1);
		int len = line.length();
		// rtrim
		while (line.charAt(len - 1) == ' ') {
			line = line.substring(0, len - 1);
			len--;
		}
		String[] nums = line.split(" ");
		double[] res = new double[nums.length];
		for (int k = 0; k < nums.length; k++) {
			if (!"".equals(nums[k])) {
				res[k] = Double.valueOf(nums[k].trim()).doubleValue();
			}
		}
		return res;
	}

	public double[] getPr() {
		return pr;
	}

	public void setPr(double[] pr) {
		this.pr = pr;
	}

	public double[][] getrV() {
		return rV;
	}

	public void setrV(double[][] rV) {
		this.rV = rV;
	}

	public double[][] getrD() {
		return rD;
	}

	public void setrD(double[][] rD) {
		this.rD = rD;
	}

	public double[][] getrV1() {
		return rV1;
	}

	public void setrV1(double[][] rV1) {
		this.rV1 = rV1;
	}

	public double[][] getDD() {
		return dD;
	}

	public void setDD(double[][] dd) {
		dD = dd;
	}

	public double[][] getDV() {
		return dV;
	}

	public void setDV(double[][] dv) {
		dV = dv;
	}

	public double[][] getDV1() {
		return dV1;
	}

	public void setDV1(double[][] dv1) {
		dV1 = dv1;
	}

	public double[][] getPd() {
		return pd;
	}

	public void setPd(double[][] pd) {
		this.pd = pd;
	}

	public double[][] getProb() {
		return prob;
	}

	public void setProb(double[][] prob) {
		this.prob = prob;
	}

	public double[] getPs() {
		return ps;
	}

	public void setPs(double[] ps) {
		this.ps = ps;
	}

	public double[][] getSD() {
		return sD;
	}

	public void setSD(double[][] sd) {
		sD = sd;
	}

	public double[][] getSV() {
		return sV;
	}

	public void setSV(double[][] sv) {
		sV = sv;
	}

	public double[][] getSV1() {
		return sV1;
	}

	public void setSV1(double[][] sv1) {
		sV1 = sv1;
	}
	
	public void setEntropyCalc(boolean value){
		entropycalc = value;
	}


	public boolean getEntropyCalc(){
		return entropycalc;
	}

}
