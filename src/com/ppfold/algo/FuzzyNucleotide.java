package com.ppfold.algo;

import java.text.DecimalFormat;
import java.util.Arrays;

public class FuzzyNucleotide {
	static DecimalFormat df = new DecimalFormat("0.00");
	
	double [] probability = new double [4]; // probability that nucleotide is an A, U, G, C
	
	public String toString()
	{
		return "("+df.format(probability[0])+", "+df.format(probability[1])+", "+df.format(probability[2])+", "+df.format(probability[3])+")";
	}
	
	public FuzzyNucleotide clone()
	{
		FuzzyNucleotide clone = new FuzzyNucleotide();
		clone.probability = Arrays.copyOf(this.probability, this.probability.length);
		return clone;		
	}
}
