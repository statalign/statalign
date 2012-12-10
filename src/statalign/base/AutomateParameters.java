package statalign.base;

import java.io.File;
import java.util.ArrayList;

import statalign.postprocess.utils.RNAFoldingTools;

/**
 * A class that contains the logic when to stop the burn-in, sampling, or 
 * how to determine the sampling rate. 
 * 
 * 
 * @author ingolfur
 *
 */


public class AutomateParameters {
	



	//private static boolean automateSamplingRate = true;
	private static boolean automateSamplingRate = false;
	private static boolean automateNumberOfSamplesToTake = false;
	//private static boolean automateBurnIn = true;
	private static boolean automateBurnIn = false;
	
	//steprate
	final static double decline = 0.0001;
	final static int right = 10;
	final static int wrong = 2;
	final static int total = right + wrong;
	final static double errorForTheAVG = 0.99;
	
	//stopsampling
	public final static double PERCENT_CONST = 0.998;
	
	//burn-in
	final static int DECLINE = 30;
	final static int NUTR_SIZE = 100;
	
	
	public static boolean shouldAutomateBurnIn() {
		return automateBurnIn;
	}
	
	public static void setAutomateBurnIn(boolean set) {
		automateBurnIn = set;
	}

	
	public static boolean shouldAutomateStepRate() {
		return automateSamplingRate;
	}

	public static void setAutomateStepRate(boolean set) {
		automateSamplingRate = set;
	}
	
	public static boolean shouldAutomateNumberOfSamples() {
		return automateNumberOfSamplesToTake;
	}

	public static void setAutomateNumberOfSamples(boolean set) {
		automateNumberOfSamplesToTake = set;
	}
	
	/**
	 * Looks either for major decline in theSpace or if we get very close to the average line. Then we break.
	 * 
	 * @param theSpace 				Gives you the average similarity between samples that have index * stepSizeOfTheBurnIn
	 * 								steps between them (theSpace[2] gives you the average similarity between samples
	 * 								that are 2 * stepSizeOfTheBurnIn steps away)
	 * @param stepSizeOfTheBurnIn   only used in the end, maybe unnecessary ?
	 * @return						the desired sampling rate (int)
	 */

	public static int getSampleRateOfTheSpace(ArrayList<Double> theSpace, int stepSizeOfTheBurnIn) {
		ArrayList<Double> der = new ArrayList<Double>(theSpace.size()-1);
		for(int i = 1; i<theSpace.size()-1; i++){
			der.add(theSpace.get(i)-theSpace.get(i+1));
		}
		int half = theSpace.size()/2;
		double avg = 0;
		int count = 0;
		for(int i =half; i<theSpace.size(); ++i ){
			avg += theSpace.get(i);
			count++;
		}
		avg = avg / count;

		

		for(int i = 1; i<(der.size()-total); ++i ){
			double f = theSpace.get(i);
			int isUprising = 0;
			for(int k = 0; k<total; ++k){
				isUprising += (der.get(i+k) < decline) ? 1 : 0;
			}
			if (isUprising >= right || f * errorForTheAVG <= avg){
				return i*stepSizeOfTheBurnIn;
			}
		}
		return 50;
	}
	
	/**
	 * If the last element in the distances variable is larger than PERCENT_CONST
	 * we stop.
	 * 
	 * @param distances 	ArrayList of distances between consecutive fuzzyAlignments
	 * @return 				True if we want to stop, else false
	 */
	
	public static boolean shouldStopSampling(ArrayList<Double> distances){
		if(distances.size() < 2){
			return false;
		}
		if(distances.get(distances.size()-1) > PERCENT_CONST){
			System.out.println("CYCLES CONVERGENCE="+distances.get(distances.size()-1) + " < " + PERCENT_CONST);
			return true;
		}
		return false;
	}

	
	/**
	 * If the logLikelihood has a major decline, we stop the burnin
	 * 
	 * @param logLikeList   List of the values where we look for the decline
	 * @return				true if we found a major decline in the logLikeList and we want to stop sampling, else false
	 */
	public static boolean shouldStopBurnIn(ArrayList<Double> logLikeList){
		if(logLikeList.size() < (NUTR_SIZE + DECLINE)){
			return false;
		}
		
		ArrayList <Double> smalllist = new ArrayList<Double>();
		ArrayList <Double> derlist = new ArrayList<Double>();
		
		for(int i =logLikeList.size()-NUTR_SIZE - DECLINE; i<logLikeList.size()-NUTR_SIZE; i++){
			double mean = 0;
			for(int j =0; j<NUTR_SIZE; j++){
				mean += logLikeList.get(i+j);
			}
		smalllist.add(mean/NUTR_SIZE);
		}

		for(int i = 0; i<smalllist.size()-1; i++){
			derlist.add(smalllist.get(i)-smalllist.get(i+1));
		}
		
		for(double i : derlist){
			if(i<0){
				return false;
			}
		}
		
		return true;
		
	}

}
