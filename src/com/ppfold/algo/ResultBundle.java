package com.ppfold.algo;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Contains final results of inside-outside algorithm. Returned to by algorithm
 * for processing. (NOTE: does not include previously excluded columns (ie.
 * those with too many gaps))
 * 
 * @author Z.Sukosd
 */
public class ResultBundle implements Serializable {

	private static final long serialVersionUID = 6577798471645335881L;
	float[][] expectationvalues;
	float[][] basepairprob;
	float[] singlebaseprob;
	char[] structure; // . = single-stranded; (,)=basepaired
	float[] reliability; // corresponding ps or pd value


	public float[][] finalmatrix;
	public double entropyVal = 0;
	public double entropyPercOfMax = 0;
	public double entropyMax = 0;
	public double reliabilityScore = 0;
	public double pairsOnlyReliabilityScore = 0;
	public double[][] phyloProbs;
	public List<Integer> leftOutColumns;

	public ResultBundle() {
	}

	public static ResultBundle tinyBundle(){
		ResultBundle tinyBundle = new ResultBundle();
		tinyBundle.structure = new char[1];
		tinyBundle.structure[0] = '.'; 
		tinyBundle.reliability = new float[1];
		tinyBundle.reliability[0] = 1;
		tinyBundle.singlebaseprob = new float[1];
		tinyBundle.singlebaseprob[0] = 1; 
		tinyBundle.basepairprob = new float[1][1]; 
		tinyBundle.basepairprob[0][0] = 0;
		tinyBundle.expectationvalues = new float[1][1]; 
		tinyBundle.expectationvalues[0][0] = 1;
		return tinyBundle;
	}
	
	public char[] getStructure() {
		return structure;
	}

	public float[][] getBasePairProb() {
		return basepairprob;
	}

	public float[][] getExpectation() {
		return expectationvalues;
	}

	public float[] getSingleBaseProb() {
		return singlebaseprob;
	}

	public float[] getReliability() {
		return reliability;
	}
	
	public float[][] getFinalMatrix()
	{
		return finalmatrix;
	}
	
	public double getPPfoldReliability()
	{
		double sum = 0;
		for(int i = 0 ; i < reliability.length ; i++)
		{
			sum += reliability[i];
		}
		return sum / ((double)reliability.length);
	}
	
	public double getPairsOnlyReliability()
	{
		double sum = 0;
		double pairs = 0;
		for(int i = 0 ; i < reliability.length ; i++)
		{
			if(structure[i] != '.')
			{
				sum += reliability[i];
				pairs++;
			}
		}
		if(pairs == 0)
		{
			return 1;
		}
		return sum / pairs;
	}
	
	public ResultBundle getSmallBundle()
	{
		ResultBundle smallBundle = new ResultBundle();
		smallBundle.structure = structure;
		smallBundle.singlebaseprob = this.singlebaseprob;
		smallBundle.entropyVal = this.entropyVal;
		smallBundle.entropyMax = this.entropyMax;
		smallBundle.entropyPercOfMax = this.entropyPercOfMax;
		smallBundle.reliabilityScore = this.getPPfoldReliability();
		smallBundle.pairsOnlyReliabilityScore = this.getPairsOnlyReliability();
		//smallBundle.phyloProbs = this.phyloProbs;
		return smallBundle;
	}
}
