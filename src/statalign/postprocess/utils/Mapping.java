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
	
    /**
     * Given an aligned sequence of length n from an alignment of length n and a n*n matrix corresponding to the alignment. This method removes gaps from the sequence to generate a sequence of length m, where m <= n and projects the specified n*n matrix
     * onto a matrix of size m*m.
     * @param alignedSequence a sequence of length n from an alignment of length n.
     * @param matrix a n*n matrix.
     * @return the projected matrix.
     */
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
}
