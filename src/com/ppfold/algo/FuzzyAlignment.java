package com.ppfold.algo;

import java.util.ArrayList;
import java.util.List;

import statalign.postprocess.utils.Mapping;

public class FuzzyAlignment {
	
	public int length;
	public List<String> sequences;
	public List<String> names;
	public List<AlignmentData> alignments;
	public List<FuzzyNucleotide[]> columns;
	public List<Integer> mapping;
	public boolean useExpectedFrequencies = true;
	
	public int getNumSequences()
	{
		return sequences.size();
	}
	
	public static FuzzyAlignment getFuzzyAlignment(List<AlignmentData> alignments)
	{
		FuzzyAlignment fuzzyAlignment = new FuzzyAlignment();
		fuzzyAlignment.sequences = alignments.get(0).sequences;
		fuzzyAlignment.names = alignments.get(0).names;
		fuzzyAlignment.length = alignments.get(0).sequences.get(0).length();
		fuzzyAlignment.alignments = alignments;
		
		int numSequences = alignments.get(0).sequences.size();
		double numAlignmentsDouble = (double) alignments.size();
		int length = alignments.get(0).sequences.get(0).length();
		fuzzyAlignment.columns = new ArrayList<FuzzyNucleotide[]>(length);
		for(int i = 0 ; i < length ; i++)
		{
			FuzzyNucleotide [] fuzzyNucleotide = new FuzzyNucleotide[numSequences];
			for(int j = 0 ; j < fuzzyNucleotide.length ; j++)
			{
				fuzzyNucleotide[j] = new FuzzyNucleotide();
			}
			fuzzyAlignment.columns.add(fuzzyNucleotide);
		}
		
		for(int al = 0 ; al < alignments.size() ; al++)
		{
			AlignmentData a = alignments.get(al);
			for(int row = 0 ; row < a.sequences.size() ; row++)
			{
				String seq = a.sequences.get(row);
				for(int col = 0 ; col < seq.length() ; col++)
				{
					double [] ambiguities = MatrixTools.createSVector(seq.charAt(col));
					for(int l = 0 ; l < ambiguities.length ; l++)
					{
						fuzzyAlignment.columns.get(col)[row].probability[l] += ambiguities[l] / numAlignmentsDouble;
					}
				}
			}
		}
		
		//System.out.println(fuzzyAlignment.length+" fuzzy\n"+fuzzyAlignment+"\n");
		return fuzzyAlignment;
	}
	
	public static double[] createSVector2(final char nt) {
		final double[] vector = MatrixTools.createVector(4, 0.01); //uncertainty in nucleotide
		//final double[] vector = createVector(4, 0.00);

		final char[] nts = MatrixTools.createNts(nt);

		double count = 0;
		for (final char n : nts) {

			if (n == 'a') {
				vector[0] = 1;
				count++;
			}
			if (n == 'u') {
				vector[1] = 1;
				count++;
			}
			if (n == 'g') {
				vector[2] = 1;
				count++;
			}
			if (n == 'c') {
				vector[3] = 1;
				count++;
			}
		}
		for(int i = 0 ; i < vector.length ; i++)
		{
			vector[i] /= count;
		}
		return vector;
	}
	
	public static FuzzyAlignment getFuzzyAlignmentAndProject(List<AlignmentData> alignments, int seqno)
	{
		ArrayList<AlignmentData> projectedAlignments = new ArrayList<AlignmentData>();
		for(int i = 0 ; i < alignments.size() ; i++)
		{
			projectedAlignments.add(projectAlignment(alignments.get(i), seqno));			
		}
		
		return getFuzzyAlignment(projectedAlignments);
	}
	
	public static FuzzyAlignment getFuzzyAlignmentAndProject(List<AlignmentData> alignments, String refSeqName)
	{
		ArrayList<AlignmentData> projectedAlignments = new ArrayList<AlignmentData>();
		for(int i = 0 ; i < alignments.size() ; i++)
		{
			projectedAlignments.add(projectAlignment(alignments.get(i), refSeqName));			
		}
		
		return getFuzzyAlignment(projectedAlignments);
	}
	
	public static AlignmentData projectAlignment(AlignmentData input, int seqno)
	{
		AlignmentData projectedAlignment = new AlignmentData();
		projectedAlignment.sequences = new ArrayList<String>();
		projectedAlignment.names = input.names;
		String refSeq = input.sequences.get(seqno);
		for(int i = 0 ; i < input.sequences.size() ; i++)
		{
			projectedAlignment.sequences.add(Mapping.projectSequence(refSeq, input.sequences.get(i), '-'));
		}
		return projectedAlignment;
	}
	
	public static AlignmentData projectAlignment(AlignmentData input, String refSeqName)
	{
		int seqno = input.names.indexOf(refSeqName);
		if(seqno == -1)
		{
			System.err.println("Could not find refseqname:" +refSeqName);
		}
		seqno = Math.max(seqno, 0);
		return projectAlignment(input, seqno);
	}
	
	public String toString()
	{
		String ret = "";
		for(int i = 0 ; i < columns.get(0).length; i++)
		{
			ret += ">"+names.get(i)+"\n"+sequences.get(i)+"\n";
			for(int j = 0 ; j < columns.size() ; j++)
			{
				ret += j+"["+columns.get(j)[i].toString()+"] ";
			}
			ret += "\n";
		}
		return ret;
	}

	public static void main(String [] args)
	{
		AlignmentData a = new AlignmentData();
		a.sequences.add("A-G");
		a.sequences.add("AAT");
		a.sequences.add("AAT");
		
		AlignmentData b = new AlignmentData();
		b.sequences.add("A-A");
		b.sequences.add("AAT");
		b.sequences.add("AAT");
		
		ArrayList<AlignmentData> alignments = new ArrayList<AlignmentData>();
		alignments.add(a);
		alignments.add(b);
		
		FuzzyAlignment fuzzy = getFuzzyAlignmentAndProject(alignments, 0);
		System.out.println(fuzzy);
	}
	
	public double [] getFrequencies(int colA, int row, boolean useMapping)
	{
		int col1 = colA;
		if(useMapping)
		{
			col1 = mapping.get(col1);
			//System.out.println(colA + " ->" + col1);
		}
		
		double [] frequencies = new double[4];
		double numAlignments = (double) alignments.size();
		
		for(int i = 0 ; i < alignments.size() ; i++)
		{
			String seq1 = alignments.get(i).sequences.get(row);
			double [] ambiguities = MatrixTools.createSVector(seq1.charAt(col1));
			for(int k = 0 ; k < ambiguities.length ; k++)
			{
				frequencies[k] += ambiguities[k]/numAlignments; 
			}
		}
		return frequencies;
	}
	
	/*public static double [][] getFrequencyPairs(List<FuzzyNucleotide[]> columns1, List<FuzzyNucleotide[]> columns2)
	{
		
	}*/
	
	
	public double [][] getFrequencyPairs(int colA, int colB, int row, boolean useMapping)
	{
		
		int col1 = colA;
		int col2 = colB;
		
		if(useMapping)
		{
			col1 = mapping.get(col1);
			col2 = mapping.get(col2);
		}
		//System.out.println("a,b="+col1+","+col2);
		double [][] nucleotidePairs = new double[4][4];
		double numAlignments = (double) alignments.size();
		
		for(int i = 0 ; i < alignments.size() ; i++)
		{
			String seq1 = alignments.get(i).sequences.get(row);
			double [] ambiguities1 = MatrixTools.createSVector(seq1.charAt(col1));
			double [] ambiguities2 = MatrixTools.createSVector(seq1.charAt(col2));
			double [][] result = new double[4][4];
			MatrixTools.multiplyVectorVector(ambiguities1, ambiguities2, result);
			
			for(int j = 0 ; j < 4 ; j++)
			{
				for(int k = 0 ; k < 4 ; k++)
				{
					nucleotidePairs[j][k] += result[j][k]/numAlignments; 
				}
			}
		}
		
		return nucleotidePairs;
	}
}


