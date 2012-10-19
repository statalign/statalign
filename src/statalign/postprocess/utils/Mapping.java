package statalign.postprocess.utils;

import java.util.Arrays;


public class Mapping {	
	/**
	 * Returns the position in the sequence without gaps.
	 * @param sequence a gapped sequence.
	 * @param gappedPos the position with gaps (just the character position).
	 * @return a position in the string ignoring gaps.
	 */
    public static int getUngappedPosition(String sequence, int gappedPos, char gapChar) {
        int ungappedPos = -1;
        int i = 0;
        int end = Math.min(sequence.length() - 1, gappedPos);
        for (i = 0; i <= end; i++) {
            if (sequence.charAt(i) != gapChar) {
                ungappedPos++;
            }
        }
        if (i == gappedPos) {
            return (ungappedPos + 1);
        }
        return ungappedPos;
    }
    
    public static float [][] projectMatrix (String alignedSequence, float [][] matrix, char gapChar)
	{		
		String ungappedSequence = alignedSequence.replaceAll("-", "");
		int projectedLength = ungappedSequence.length();
		int [] gappedToUngapped = new int[alignedSequence.length()];
		int [] ungappedToGapped = new int[projectedLength];
		Arrays.fill(gappedToUngapped, -1);
		Arrays.fill(ungappedToGapped, -1);
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{			
			if(gappedToUngapped[i] == -1)
			{
				int x = getUngappedPosition(alignedSequence, i, '-');
				gappedToUngapped[i] = x;
			}
		}
		
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
			{
				ungappedToGapped[gappedToUngapped[i]] = i;
			}
		}
		
		/*
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			System.out.println(i + " -> " + gappedToUngapped[i]);
		}
		System.out.println();
		for(int i = 0 ; i < ungappedToGapped.length ; i++)
		{	
			System.out.println(i + " -> " + ungappedToGapped[i]);
		}*/
		
		float [][] projectedMatrix = new float[projectedLength][projectedLength];
		for(int i = 0 ; i < matrix.length ; i++)
		{
			for(int j = 0 ; j < matrix[0].length ; j++)
			{
				if(gappedToUngapped[i] != -1 && gappedToUngapped[j] != -1)
				{
					projectedMatrix[gappedToUngapped[i]][gappedToUngapped[j]] = matrix[i][j];
				}
			}
		}
		
		/*
		
		// renormalize
		double matrixSum = 0;
		double projectSum = 0;
		for(int i = 0 ; i < matrix.length ; i++)
		{
			for(int j = 0 ; j < matrix[0].length ; j++)
			{
				matrixSum += matrix[i][j];
			}
		}
		
		for(int i = 0 ; i < projectedMatrix.length ; i++)
		{
			for(int j = 0 ; j < projectedMatrix[0].length ; j++)
			{
				projectSum += projectedMatrix[i][j];
			}
		}
		double scaleFactor = matrixSum / projectSum;
		for(int i = 0 ; i < projectedMatrix.length ; i++)
		{
			for(int j = 0 ; j < projectedMatrix[0].length ; j++)
			{
				projectedMatrix[i][j] = (float)(scaleFactor*projectedMatrix[i][j]);
				if(projectedMatrix[i][j] > 1)
				{
					System.out.println("Projected matrix > 1 " + i + ", " + j + " : " + projectedMatrix[i][j]);
				}
			}
		}
		*/
		
		return projectedMatrix;
	}
    
    public static int [] getProjectionIndices(String alignedSequence, char gapChar)
    {
    	String ungappedSequence = alignedSequence.replaceAll("-", "");
   		int projectedLength = ungappedSequence.length();
   		int [] gappedToUngapped = new int[alignedSequence.length()];
   		int [] ungappedToGapped = new int[projectedLength];
   		Arrays.fill(gappedToUngapped, -1);
   		Arrays.fill(ungappedToGapped, -1);
   		for(int i = 0 ; i < gappedToUngapped.length ; i++)
   		{			
   			if(gappedToUngapped[i] == -1)
   			{
   				int x = getUngappedPosition(alignedSequence, i, '-');
   				gappedToUngapped[i] = x;
   			}
   		}
   		
   		for(int i = 0 ; i < gappedToUngapped.length ; i++)
   		{	
   			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
   			{
   				ungappedToGapped[gappedToUngapped[i]] = i;
   			}
   		}
   		
   		return gappedToUngapped;
    }
    
    public static double [][] projectMatrix (String alignedSequence, double [][] matrix, char gapChar)
   	{		
   		String ungappedSequence = alignedSequence.replaceAll("-", "");
   		int projectedLength = ungappedSequence.length();
   		int [] gappedToUngapped = new int[alignedSequence.length()];
   		int [] ungappedToGapped = new int[projectedLength];
   		Arrays.fill(gappedToUngapped, -1);
   		Arrays.fill(ungappedToGapped, -1);
   		for(int i = 0 ; i < gappedToUngapped.length ; i++)
   		{			
   			if(gappedToUngapped[i] == -1)
   			{
   				int x = getUngappedPosition(alignedSequence, i, '-');
   				gappedToUngapped[i] = x;
   			}
   		}
   		
   		for(int i = 0 ; i < gappedToUngapped.length ; i++)
   		{	
   			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
   			{
   				ungappedToGapped[gappedToUngapped[i]] = i;
   			}
   		}
   		
   		/*
   		for(int i = 0 ; i < gappedToUngapped.length ; i++)
   		{	
   			System.out.println(i + " -> " + gappedToUngapped[i]);
   		}
   		System.out.println();
   		for(int i = 0 ; i < ungappedToGapped.length ; i++)
   		{	
   			System.out.println(i + " -> " + ungappedToGapped[i]);
   		}*/
   		
   		double [][] projectedMatrix = new double[projectedLength][projectedLength];
   		for(int i = 0 ; i < matrix.length ; i++)
   		{
   			for(int j = 0 ; j < matrix[0].length ; j++)
   			{
   				if(gappedToUngapped[i] != -1 && gappedToUngapped[j] != -1)
   				{
   					projectedMatrix[gappedToUngapped[i]][gappedToUngapped[j]] = matrix[i][j];
   				}
   			}
   		}
   		
   		/*
   		
   		// renormalize
   		double matrixSum = 0;
   		double projectSum = 0;
   		for(int i = 0 ; i < matrix.length ; i++)
   		{
   			for(int j = 0 ; j < matrix[0].length ; j++)
   			{
   				matrixSum += matrix[i][j];
   			}
   		}
   		
   		for(int i = 0 ; i < projectedMatrix.length ; i++)
   		{
   			for(int j = 0 ; j < projectedMatrix[0].length ; j++)
   			{
   				projectSum += projectedMatrix[i][j];
   			}
   		}
   		double scaleFactor = matrixSum / projectSum;
   		for(int i = 0 ; i < projectedMatrix.length ; i++)
   		{
   			for(int j = 0 ; j < projectedMatrix[0].length ; j++)
   			{
   				projectedMatrix[i][j] = (float)(scaleFactor*projectedMatrix[i][j]);
   				if(projectedMatrix[i][j] > 1)
   				{
   					System.out.println("Projected matrix > 1 " + i + ", " + j + " : " + projectedMatrix[i][j]);
   				}
   			}
   		}
   		*/
   		
   		return projectedMatrix;
   	}
	
    /**
     * Given an aligned sequence of length n from an alignment of length n and a n*n matrix corresponding to the alignment. This method removes gaps from the sequence to generate a sequence of length m, where m <= n and projects the specified n*n matrix
     * onto a matrix of size m*m.
     * @param alignedSequence a sequence of length n from an alignment of length n.
     * @param matrix a n*n matrix.
     * @return the projected matrix.
     */
	public static float [][] projectMatrix2 (String alignedSequence, float [][] matrix, char gapChar)
	{		
		String ungappedSequence = alignedSequence.replaceAll("-", "");
		int projectedLength = ungappedSequence.length();
		int [] gappedToUngapped = new int[alignedSequence.length()];
		int [] ungappedToGapped = new int[projectedLength];
		Arrays.fill(gappedToUngapped, -1);
		Arrays.fill(ungappedToGapped, -1);
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{			
			if(gappedToUngapped[i] == -1)
			{
				int x = getUngappedPosition(alignedSequence, i, '-');
				gappedToUngapped[i] = x;
			}
		}
		
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
			{
				ungappedToGapped[gappedToUngapped[i]] = i;
			}
		}
		
		/*
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			System.out.println(i + " -> " + gappedToUngapped[i]);
		}
		System.out.println();
		for(int i = 0 ; i < ungappedToGapped.length ; i++)
		{	
			System.out.println(i + " -> " + ungappedToGapped[i]);
		}*/
		
		float [][] projectedMatrix = new float[projectedLength][projectedLength];
		for(int i = 0 ; i < projectedMatrix.length ; i++)
		{
			for(int j = 0 ; j < projectedMatrix[0].length ; j++)
			{
				if(ungappedToGapped[i] != -1 && ungappedToGapped[j] != -1)
				{
					projectedMatrix[i][j] = matrix[ungappedToGapped[i]][ungappedToGapped[j]];
				}
			}
		}
		
		return projectedMatrix;
	}
	
	 public static float [] projectarray (String alignedSequence, float [] array, char gapChar)
	 {		
		String ungappedSequence = alignedSequence.replaceAll("-", "");
		int projectedLength = ungappedSequence.length();
		int [] gappedToUngapped = new int[alignedSequence.length()];
		int [] ungappedToGapped = new int[projectedLength];
		Arrays.fill(gappedToUngapped, -1);
		Arrays.fill(ungappedToGapped, -1);
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{			
			if(gappedToUngapped[i] == -1)
			{
				int x = getUngappedPosition(alignedSequence, i, '-');
				gappedToUngapped[i] = x;
			}
		}
		
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
			{
				ungappedToGapped[gappedToUngapped[i]] = i;
			}
		}
		
		/*
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			System.out.println(i + " -> " + gappedToUngapped[i]);
		}
		System.out.println();
		for(int i = 0 ; i < ungappedToGapped.length ; i++)
		{	
			System.out.println(i + " -> " + ungappedToGapped[i]);
		}*/
		
		float [] projectedArray = new float[projectedLength];
		for(int i = 0 ; i < array.length ; i++)
		{
			if(gappedToUngapped[i] != -1)
			{
				projectedArray[gappedToUngapped[i]] = array[i];
			}
		}
		
		return projectedArray;
	}
	 
	 public static String projectSequence2 (String alignedSequence, String sequenceFromAlignment, char gapChar)
	 {		
		String ungappedSequence = alignedSequence.replaceAll(gapChar+"", "");
		int projectedLength = ungappedSequence.length();
		/*int [] gappedToUngapped = new int[alignedSequence.length()];
		int [] ungappedToGapped = new int[projectedLength];
		Arrays.fill(gappedToUngapped, -1);
		Arrays.fill(ungappedToGapped, -1);
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{			
			if(gappedToUngapped[i] == -1)
			{
				int x = getUngappedPosition(alignedSequence, i, gapChar);
				gappedToUngapped[i] = x;
			}
		}
		
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
			{
				ungappedToGapped[gappedToUngapped[i]] = i;
			}
		}
		*/
		/*
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			System.out.println(i + " -> " + gappedToUngapped[i]);
		}
		System.out.println();
		for(int i = 0 ; i < ungappedToGapped.length ; i++)
		{	
			System.out.println(i + " -> " + ungappedToGapped[i]);
		}*/
		String gaps = "";
		for(int i = 0 ; i < projectedLength ; i++)
		{
			gaps += "-";
		}
		
		StringBuffer projectedSequence = new StringBuffer(gaps);
		int k = 0;
		for(int i = 0 ; i < sequenceFromAlignment.length() ; i++)
		{
			if(alignedSequence.charAt(i) != gapChar)
			{
				projectedSequence.setCharAt(k, sequenceFromAlignment.charAt(i));
				k++;
			}
			
		}
		
		return projectedSequence.toString();
	}
	 
	 public static String projectSequence (String alignedSequence, String sequenceFromAlignment, char gapChar)
	 {		
		String ungappedSequence = alignedSequence.replaceAll(gapChar+"", "");
		int projectedLength = ungappedSequence.length();
		int [] gappedToUngapped = new int[alignedSequence.length()];
		int [] ungappedToGapped = new int[projectedLength];
		Arrays.fill(gappedToUngapped, -1);
		Arrays.fill(ungappedToGapped, -1);
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{			
			if(gappedToUngapped[i] == -1)
			{
				int x = getUngappedPosition(alignedSequence, i, gapChar);
				gappedToUngapped[i] = x;
			}
		}
		
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
			{
				ungappedToGapped[gappedToUngapped[i]] = i;
			}
		}
		
		/*
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			System.out.println(i + " -> " + gappedToUngapped[i]);
		}
		System.out.println();
		for(int i = 0 ; i < ungappedToGapped.length ; i++)
		{	
			System.out.println(i + " -> " + ungappedToGapped[i]);
		}*/
		String gaps = "";
		for(int i = 0 ; i < projectedLength ; i++)
		{
			gaps += "-";
		}
		
		StringBuffer projectedSequence = new StringBuffer(gaps);
		int k = 0;
		for(int i = 0 ; i < sequenceFromAlignment.length() ; i++)
		{
			if(alignedSequence.charAt(i) != gapChar)
			{
				//projectedArray[gappedToUngapped[i]] = array[i];
				projectedSequence.setCharAt(gappedToUngapped[i], sequenceFromAlignment.charAt(i));
				//k++;
			}
			
		}
		
		return projectedSequence.toString();
	}
	
	public static float [] projectArray2 (String alignedSequence, float [] array, char gapChar)
	{		
		String ungappedSequence = alignedSequence.replaceAll("-", "");
		int projectedLength = ungappedSequence.length();
		int [] gappedToUngapped = new int[alignedSequence.length()];
		int [] ungappedToGapped = new int[projectedLength];
		Arrays.fill(gappedToUngapped, -1);
		Arrays.fill(ungappedToGapped, -1);
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{			
			if(gappedToUngapped[i] == -1)
			{
				int x = getUngappedPosition(alignedSequence, i, '-');
				gappedToUngapped[i] = x;
			}
		}
		
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
			{
				ungappedToGapped[gappedToUngapped[i]] = i;
			}
		}
		
		/*
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			System.out.println(i + " -> " + gappedToUngapped[i]);
		}
		System.out.println();
		for(int i = 0 ; i < ungappedToGapped.length ; i++)
		{	
			System.out.println(i + " -> " + ungappedToGapped[i]);
		}*/
		
		float [] projectedArray = new float[projectedLength];
		for(int i = 0 ; i < projectedArray.length ; i++)
		{
				if(ungappedToGapped[i] != -1 && ungappedToGapped[i] < array.length)
				{
					projectedArray[i] = array[ungappedToGapped[i]];
				}
		}
		
		return projectedArray;
	}
	

	
	public static int [] getUngappedToGappedMapping(String alignedSequence)
	{
		String ungappedSequence = alignedSequence.replaceAll("-", "");
		int projectedLength = ungappedSequence.length();
		int [] gappedToUngapped = new int[alignedSequence.length()];
		int [] ungappedToGapped = new int[projectedLength];
		Arrays.fill(gappedToUngapped, -1);
		Arrays.fill(ungappedToGapped, -1);
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{			
			if(gappedToUngapped[i] == -1)
			{
				int x = getUngappedPosition(alignedSequence, i, '-');
				gappedToUngapped[i] = x;
			}
		}
		
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
			{
				ungappedToGapped[gappedToUngapped[i]] = i;
			}
		}
		
		
		return ungappedToGapped;
	}
	
	public static int [] getGappedToUngappedMapping(String alignedSequence)
	{
		String ungappedSequence = alignedSequence.replaceAll("-", "");
		int projectedLength = ungappedSequence.length();
		int [] gappedToUngapped = new int[alignedSequence.length()];
		int [] ungappedToGapped = new int[projectedLength];
		Arrays.fill(gappedToUngapped, -1);
		Arrays.fill(ungappedToGapped, -1);
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{			
			if(gappedToUngapped[i] == -1)
			{
				int x = getUngappedPosition(alignedSequence, i, '-');
				gappedToUngapped[i] = x;
			}
		}
		
		for(int i = 0 ; i < gappedToUngapped.length ; i++)
		{	
			if(gappedToUngapped[i] != -1 && ungappedToGapped[gappedToUngapped[i]] == -1)
			{
				ungappedToGapped[gappedToUngapped[i]] = i;
			}
		}
		
		
		return gappedToUngapped;
	}
	
	public static int [] projectPairedSites(String alignedSequence, int [] pairedSites)
	{
		int [] ungappedToGapped = Mapping.getUngappedToGappedMapping(alignedSequence);
		int [] gappedToUngapped = Mapping.getGappedToUngappedMapping(alignedSequence);
		
		int [] projectedPairedSites = new int[ungappedToGapped.length];
		for(int i = 0 ; i < pairedSites.length ; i++)
		{
			if(pairedSites[i] != 0) // if paired, map
			{
				int x = gappedToUngapped[i];
				if(x != -1);
				{
					int y = Math.max(0, gappedToUngapped[pairedSites[i]-1]) + 1;
					projectedPairedSites[x] = y;
				}
			}
		}
		return projectedPairedSites;
	}
}
