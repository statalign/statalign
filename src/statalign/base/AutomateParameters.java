package statalign.base;

import java.util.ArrayList;

public class AutomateParameters {
	
	

	public static void main(String[] args){
		double[][] a = new double[2][2];
		double[][] b = new double[2][2];

		a[0][0] = 0;
		a[0][1] = 1;
		a[1][0] = 2;
		a[1][1] = 3;

		b[0][0] = 0.1;
		b[0][1] = 0.2;
		b[1][0] = 0.3;
		b[1][1] = 0.4;

		//System.out.println(hasConverged(a,b));


	}



	private static boolean automateSamplingRate = true;
	private static boolean automateNumberOfSamplesToTake = true;
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
		System.out.println(theSpace);
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

		double decline = 0.0001;
		int right = 10;
		int wrong = 2;
		int total = right + wrong;

		for(int i = 1; i<(der.size()-total); ++i ){
			double f = theSpace.get(i);
			int isUprising = 0;
			for(int k = 0; k<total; ++k){
				isUprising += (der.get(i+k) < decline) ? 1 : 0;
			}
			if (isUprising >= right || f * 0.99 <= avg){
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
		return logLikeList.size() > 10000;
	}

}
