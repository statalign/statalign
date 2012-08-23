package statalign.distance;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import statalign.postprocess.utils.RNAFoldingTools;

/**
 * 
 * A class that calculates the distance and similarity score between two alignments
 * 
 * @author Ingolfur
 *
 */

public class Distance {


	/**
	 * Given a set of sequences and an AMA score, returns the corresponding distance.
	 * @param sequences  sequences (possibly having gaps)
	 * @param amaScore	 similarity score between 0 and 1
	 * @return  		 the distance between the alignments
	 */
	public static double amaScoreToDistance(ArrayList<String> sequences, double amaScore)
	{
		double lengthSum = 0;
		for(int i = 0 ; i < sequences.size() ; i++)
		{
			lengthSum += sequences.get(i).replaceAll("-", "").length();
		}
		double ksub1 = sequences.size() - 1;

		return (1-amaScore)*ksub1*lengthSum;
	}
	/**
	 * So we can use Strings instead of list
	 * 
	 * @see #distance(ArrayList, ArrayList)
	 * 
	 */

	public static int distance(String [] A, String [] B){
		ArrayList<String> arListA = new ArrayList <String> (Arrays.asList(A));
		ArrayList<String> arListB = new ArrayList <String> (Arrays.asList(B));
		return Distance.distance(arListA, arListB);

	}
	
	/**
	 * Calculates the distance between two alignments, if the function
	 * returns 0, the alignments are the same. Higher number means the
	 * alignments are more distant.
	 * 
	 * 
	 * @param A     The first alignment
	 * @param B		The second alignment	
	 * @return		Integer between 0 and infinity
	 */

	public static int distance(ArrayList<String>  A, ArrayList<String>  B){
		ArrayList<String> cloneA = new ArrayList<String>();
		ArrayList<String> cloneB = new ArrayList<String>();

		for(String i : A){
			cloneA.add(i);
		}

		for(String i : B){
			cloneB.add(i);
		}

		sortSeq(cloneA,cloneB);
		String [] tempA = new String[2]; 
		String [] tempB = new String[2];
		int d = 0;
		int k = A.size(); //This says how many sequences the alignments have
		for(int i=0; i<k-1; ++i){
			for(int j=i+1; j<k; ++j){
				tempA[0] = A.get(i);
				tempA[1] = A.get(j);
				tempB[0] = B.get(i);
				tempB[1] = B.get(j);
				d += Distance.distanceBetweenTwoAlignments(tempA, tempB);
			}
		}

		return d;
	}

	private static void sortSeq(ArrayList<String>  A, ArrayList<String>  B){

		List< Pair < String  , Integer> > seq1 = new ArrayList< Pair < String , Integer> >();
		List< Pair < String  , Integer> > seq2 = new ArrayList< Pair < String  , Integer> >();
		ArrayList<String> newA = new ArrayList<String>();
		ArrayList<String> newB = new ArrayList<String>();
		for(int i = 0; i<A.size(); ++i){
			seq1.add( new Pair<String,Integer >  (A.get(i).replaceAll("-",""),i));
			seq2.add( new Pair<String,Integer >  (B.get(i).replaceAll("-",""),i));
		}
		Collections.sort(seq1, new PairCompare());
		Collections.sort(seq2, new PairCompare());

		for(int i = 0; i<A.size(); ++i){
			int index1 =seq1.get(i).getRight();
			int index2 =seq2.get(i).getRight();

			newA.add(A.get(index1));
			newB.add(B.get(index2));
		}
		for(int i = 0; i<A.size(); ++i){
			A.set(i, newA.get(i));
			B.set(i, newB.get(i));
		}
	}

	private static int distanceBetweenTwoAlignments(String[] A, String[] B){
		String seqA1 = A[0].replaceAll("-", "");
		String seqA2 = A[1].replaceAll("-", "");

		int homoSize = Distance.HomoUnionSize(A, B);
		int delandInsertSize = deletionUnionAndInsertionUnionSize(A,B); 

		int dist = seqA1.length() + seqA2.length() - 2 * homoSize - delandInsertSize;
		return dist;

	}

	private static int HomoUnionSize(String[] A, String[] B) {
		ArrayList<Pair<Integer,Integer>> HomoForA = new ArrayList<Pair<Integer,Integer>>();
		int seqIndex1 = 0;
		int seqIndex2 = 0;
		char seq1Char;
		char seq2Char;
		for(int i =0; i<A[0].length(); ++i ){
			seq1Char = A[0].charAt(i);
			seq2Char = A[1].charAt(i);

			if(seq2Char != '-' && seq1Char != '-') {
				HomoForA.add( new Pair<Integer,Integer>(seqIndex1,seqIndex2));
			}

			if(seq1Char != '-'){
				seqIndex1++;
			}


			if(seq2Char != '-'){
				seqIndex2++;
			}

		}

		ArrayList<Pair<Integer,Integer>> HomoForB = new ArrayList<Pair<Integer,Integer>>();
		seqIndex1 = 0;
		seqIndex2 = 0;
		for(int i =0; i<B[0].length(); ++i ){
			seq1Char = B[0].charAt(i);
			seq2Char = B[1].charAt(i);

			if(seq2Char != '-' && seq1Char != '-') {
				HomoForB.add( new Pair<Integer,Integer>(seqIndex1,seqIndex2));
			}

			if(seq1Char != '-'){
				seqIndex1++;
			}

			if(seq2Char != '-'){
				seqIndex2++;
			}
		}

		ArrayList<Pair<Integer,Integer>> HomoUnion = new ArrayList<Pair<Integer,Integer>>();
		for(int j=0;j<HomoForA.size();++j)
		{	
			Pair<Integer,Integer> common = HomoForA.get(j);
			if(HomoForB.contains(common)){
				HomoUnion.add(common);
			}
		}
		return HomoUnion.size();

	}

	private static int deletionUnionAndInsertionUnionSize(String[] A, String[] B) {
		ArrayList<Integer> deletionForA = new ArrayList<Integer>();
		ArrayList<Integer> insertionForA = new ArrayList<Integer>();
		char seq1Char;
		char seq2Char;
		int seqIndex1 = 0;
		int seqIndex2 = 0;
		for(int i =0; i<A[0].length(); ++i ){
			seq1Char = A[0].charAt(i);
			seq2Char = A[1].charAt(i);

			if(seq1Char == '-' && seq2Char == '-'){
				continue;
			}

			if(seq1Char == '-'){
				insertionForA.add(seqIndex2);
			}
			else{
				seqIndex1++;
			}

			if(seq2Char == '-'){
				deletionForA.add(seqIndex1-1);
			}
			else{
				seqIndex2++;
			}


		}

		seqIndex1=0;
		seqIndex2=0;
		ArrayList<Integer> deletionForB = new ArrayList<Integer>();
		ArrayList<Integer> insertionForB = new ArrayList<Integer>();
		for(int i =0; i<B[0].length(); ++i ){
			seq1Char = B[0].charAt(i);
			seq2Char = B[1].charAt(i);

			if(seq1Char == '-' && seq2Char == '-'){
				continue;
			}


			if(seq2Char == '-'){
				deletionForB.add(seqIndex1);
			}
			else{
				seqIndex2++;
			}

			if(seq1Char == '-'){
				insertionForB.add(seqIndex2-1);
			}
			else{
				seqIndex1++;
			}
		}



		ArrayList<Integer> deletionUnion = new ArrayList<Integer>();
		for(int j=0;j<deletionForA.size();++j)
		{	
			int common = deletionForA.get(j);
			if(deletionForB.contains(common)){
				deletionUnion.add(common);
			}
		}

		ArrayList<Integer> insertionUnion = new ArrayList<Integer>();
		for(int j=0;j<insertionForA.size();++j)
		{	
			int common = insertionForA.get(j);
			if(insertionForB.contains(common)){
				insertionUnion.add(common);
			}
		}
		return deletionUnion.size() + insertionUnion.size();
	}
	
	/**
	 * 
	 * The same as {@link #AMA(ArrayList, ArrayList)} but here you can use arrays instead of lists
	 * 
	 */

	public static double AMA(String[] A, String[] B){
		ArrayList<String> arListA = new ArrayList <String> (Arrays.asList(A));
		ArrayList<String> arListB = new ArrayList <String> (Arrays.asList(B));
		return Distance.AMA(arListA, arListB);
	}
	
	
	/**
	 * Does the same as the {@link Distance#distance(ArrayList, ArrayList)} but outputs a more
	 * meaningful number since it is between 0 and 1. 1 means we have the same alignments while 0 means they are as distance as possible.
	 * 
	 * @param A		The first alignment
	 * @param B		The second alignment
	 * @return		Number between 0 and 1
	 */
	public static double AMA(ArrayList<String>  A, ArrayList<String>  B){
		int sum = 0;
		for(int i = 0; i<A.size(); ++i){
			sum += A.get(i).replaceAll("-","").length();
		}
		return 1 - Distance.distance(A, B)/ (double)( (A.size()-1) *sum );
	}
	
	/**
	 * Overloading the sequenceSimilarityScore method so one can use arrays instead of lists
	 */

	public static double sequenceSimilarityScore(String []  A){
		ArrayList<String> temp = new ArrayList<String>(Arrays.asList(A)); 
		return sequenceSimilarityScore(temp);
	}
	
	/**
	 * Calculates the similarity of the sequences. 
	 * 
	 * @param A	List of strings where each string is a sequence ()
	 * @return 	A double between 0 and 1, 1 means all sequences are the same. 0 they are as different as possible 
	 */
	public static double sequenceSimilarityScore(List<String>  A){
		double score = 0;
		int count = 0;
		int tempScore;
		int length;
		for(int i =0; i<A.size()-1; ++i){
			for(int j =i+1; j<A.size(); ++j){
				String seq1 = A.get(i);
				String seq2 = A.get(j);
				tempScore = 0;
				length = 0;
				count++;
				for(int k = 0; k<seq1.length(); ++k){
					if (seq1.charAt(k) == '.' || seq2.charAt(k) == '.'){
						continue;
					}
					else if(seq1.charAt(k) == seq2.charAt(k)){
						tempScore++;
						length++;
					}else{
						length++;
					}
				}
				score += tempScore / (double)length;
			}
		}
		return score /count;
	}
	
	/**
	 * Used when determining the step rate in the automation of the MCMC. Gives you
	 * a much better idea how the MCMC progress behaves that is, how the space looks like
	 * on average.
	 * 
	 * 
	 * @param allAlignments	List of alignments
	 * @return	A list of average similarity between all possible pairs of alignments
	 */

	public static ArrayList<Double> spaceAMA(ArrayList<String[]> allAlignments){
		ArrayList<Double> sum = new ArrayList<Double>(allAlignments.size()+1);
		//Initialising sum
		for(int i =0; i<allAlignments.size(); ++i){
			sum.add(0.0);
		}

		for(int i = 0; i<allAlignments.size()-1; ++i){
			for(int j = i+1; j<allAlignments.size(); ++j){
				String[] first = allAlignments.get(i);
				String[] second = allAlignments.get(j);
				sum.set(j-i, sum.get(j-i)+Distance.AMA(first, second));
			}
		} 

		for(int i = 0; i<allAlignments.size(); ++i){
			sum.set(i, sum.get(i)/(allAlignments.size() - i  ) );
		}
		return sum;
	}
}

