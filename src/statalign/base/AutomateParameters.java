package statalign.base;

import java.io.File;
import java.util.ArrayList;

import statalign.postprocess.utils.RNAFoldingTools;

public class AutomateParameters {
	



	private static boolean automateSamplingRate = true;
	private static boolean automateNumberOfSamplesToTake = false;
	private static boolean automateBurnIn = true;
	

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

	public static int getSampleRateOfTheSpace(ArrayList<Double> theSpace, int stepSizeOfTheBurnIn2) {
		final double decline = 0.0001;
		final int right = 10;
		final int wrong = 2;
		final int total = right + wrong;
		final double error = 0.99;
		
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
			if (isUprising >= right || f * error <= avg){
				return i*stepSizeOfTheBurnIn2;
			}
		}
		return 1;
	}
	
	public static boolean shouldStopSampling(ArrayList<Double> distances){
		double PERCENT_CONST = 0.999;
		if(distances.size() < 2){
			return false;
		}
		if(distances.get(distances.size()-1) > PERCENT_CONST){
			return true;
		}
		return false;
	}

	
	public static boolean shouldStopBurnIn(ArrayList<Double> logLikeList){
		final int DECLINE = 30;
		final int NUTR_SIZE = 100;
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
