package com.ppfold.algo.extradata;


import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;
import java.math.*;

import com.ppfold.algo.MatrixTools;




/**
 * Static class to define the model for how to deal with SHAPE data as well as container for the actual data 
 * 
 * @author Z.Sukosd
 */


public class SHAPEDataFctDiffBp implements ExtraData {
	
	//the percent error assumed in each individual SHAPE measurement
	//(this is used to define the interval to calculate probability integrals) 
	
	private int type = 0;
	
	private float[] data; //contains the actual measurement numbers (data)
	
	//Probabilities must always be given as P(data|unpaired) or P(data|paired).  
	private float[] dataProbGivenInnerPair; //P(data|model), inner pair case
	private float[] dataProbGivenOuterPair; //P(data|model), outer pair case
	private float[] dataProbGivenUnpaired; //P(data|model), unpaired case
	
	//Unpaired parameter fit
	private final float lamb = 0.681211f;
	
	private final float k_stacked =  0.852665f;
	private final float sigma_stacked = 0.05284387f;
	private final float mu_stacked = 0.0451519f;
	
	private final float k_helixend =  0.737508f;
	private final float sigma_helixend = 0.125595f;
	private final float mu_helixend = 0.104181f;
	
//	
//	private final float k_stacked =  0.815588f;
//	private final float sigma_stacked = 0.0526406f;
//	private final float mu_stacked = 0.041221f;
//	
//	private final float k_helixend =  0.791414f;
//	private final float sigma_helixend = 0.121313f;
//	private final float mu_helixend = 0.0983816f;
//	
//	private final float k_mixedpair =  0.895341f;
//	private final float sigma_mixedpair = 0.0712473f;
//	private final float mu_mixedpair = 0.054239f;

	public SHAPEDataFctDiffBp(){
	}
	
	public int getType(){
		return type; 
	}

	public boolean isEmpty(int i)
	{
		return data[i]==-999;
	}
	public float getProbabilityGivenOuterPaired(int position1, int position2) {
		return getProbabilityGivenOuterPair(position1)*getProbabilityGivenOuterPair(position2);
	}
	public float getProbabilityGivenInnerPaired(int position1, int position2) {
		return getProbabilityGivenInnerPair(position1)*getProbabilityGivenInnerPair(position2);
	}
	
	public float getProbabilityGivenInnerPair(int n){
		return dataProbGivenInnerPair[n];  
	}
	public float getProbabilityGivenOuterPair(int n){
		return dataProbGivenOuterPair[n];
	}
	public float getProbabilityGivenUnpaired(int n){
		return dataProbGivenUnpaired[n];
	}
	
	private float expCDF(float x, float lamb){
		return 1f - (float)Math.exp(-x/lamb);
	}
	
	private float gevCDF(float x, float xi, float sigma, float mu){
		float oneoverxi = 1f/xi;
		return (float)Math.exp(-1f*Math.pow(1 + xi*(x-mu)/sigma,-oneoverxi) ); 		
	}
	
	private float expPDF(float x, float lamb){
		return (1/lamb) * (float)Math.exp(-x/lamb);
	}
	
	private float gevPDF(float x, float xi, float sigma, float mu){
		float oneoverxi = 1f/xi;
		float oneoversigma = 1f/sigma;
		float xminusmu = x-mu;
		return oneoversigma*(float)Math.exp(-1*Math.pow(1+xi*xminusmu*oneoversigma,-1*oneoverxi))*
			(float)Math.pow((1+xi*xminusmu*oneoversigma),-1-oneoverxi); 		
	}
	
	private float calcUnpaired(float f) {
		//Calculates the numbers on the basis of distribution fits
		return expPDF(f,lamb)/10f;
	}

	private float calcOuterPair(float f) {
		// Calculates the numbers on the basis of distribution fits
		return gevPDF(f,k_helixend, sigma_helixend, mu_helixend)/10f;
	}

	private float calcInnerPair(float f) {
		// Calculates the numbers on the basis of distribution fits
		return gevPDF(f,k_stacked, sigma_stacked, mu_stacked)/10f;
		//return gevPDF(f,k_mixedpair, sigma_mixedpair, mu_mixedpair);
	}
	

	public void importData(String filename, int sequencelength) throws Exception{
		//read the SHAPE data 
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			System.err.println("SHAPE sequence input file " + filename + " could not be read!");
			throw new IOException(e);
		}
		readData_toStream(stream, sequencelength);	
	}
	
	public void readData_toStream(BufferedInputStream stream, int sequencelength) throws Exception{
		//the CLC import plugin uses the following line:  
		//public ClcObject[] doImport(BufferedInputStream stream, String name, Activity activity) throws IOException, ParseException, PersistenceException {
		//the same stream will be passed here to ease the import process later.
		
		//create the data containers
		this.data = new float[sequencelength];
		this.dataProbGivenInnerPair = new float[sequencelength];
		this.dataProbGivenOuterPair = new float[sequencelength];
		this.dataProbGivenUnpaired = new float[sequencelength];
		if(stream!=null){
			String line = "";
			try{
				int l = stream.available();
				byte[] bytes = new byte[l];
				stream.read(bytes);
				stream.close();
				String data_string = new String(bytes);
				String[] lines = data_string.split("\n");
				int data_size = lines.length;
				int [] readdata_index = new int[data_size];
				float [] readdata_data = new float[data_size];
				// Create a pattern to match different kinds of separators
		        Pattern p = Pattern.compile("[,\\s]+");
		        //DO NOT ignore first line... 
				for(int i = 0; i<lines.length; i++){
					line = lines[i];
					String splitline[] = p.split(line.trim());
					if(splitline.length == 2){
						readdata_index[i] = Integer.valueOf(splitline[0])-1; //SHAPE data files are numbered from 1; Java numbers from 0
						readdata_data[i] = Float.valueOf(splitline[1]);
					}
				}

				//initialize all data values to -999
				for(int i = 0; i<this.data.length; i++){
					data[i] = -999;
				}
				//set the appropriate data values to the read ones
				for(int i = 0; i<data_size; i++){
					//
					data[readdata_index[i]] = readdata_data[i];
				}
				
				//calculate probabilities 
				calcProbabilities();
				
				//normalize probabilities
				//normalizeProbabilities(); 
				
			}
			catch(Exception e){
				System.err.println("An exception occured while attempting to read or interpret the SHAPE sequence data. ");
				throw new Exception(e);
			}
		}
		else{
			System.err.println("Input stream was null, SHAPE sequence data could not be loaded.");
		}
	}
	
	private void calcProbabilities(){
		System.out.println("Calculating probs");
		int n = data.length;
		for(int i = 0; i<n; i++){
			//For each data point, calculate the probability of data given paired/unpaired. 
			//Assume a proportional error of ERROR_FRACTION in each measurement.
			if(data[i] == -999){
				//ignore 
				dataProbGivenInnerPair[i] = 1;
				dataProbGivenOuterPair[i] = 1;
				dataProbGivenUnpaired[i] = 1;	
				continue; 
			}

			dataProbGivenInnerPair[i] = calcInnerPair(data[i]);
			dataProbGivenOuterPair[i] = calcOuterPair(data[i]);
			dataProbGivenUnpaired[i] = calcUnpaired(data[i]);
		}
	}
	
	private void normalizeProbabilities(){
		System.out.println("Normalizing probs");
		int n = data.length;
		for(int i = 0; i<n; i++){
			//For each data point, normalize the three probabilities 
			//so they are all under 1. 
			//This is so we don't get overflow errors
			//(it's the ratios that matter anyway)
			float sum = dataProbGivenInnerPair[i] + dataProbGivenOuterPair[i] + dataProbGivenUnpaired[i];
			dataProbGivenInnerPair[i] = dataProbGivenInnerPair[i]/sum;
			dataProbGivenOuterPair[i] = dataProbGivenOuterPair[i]/sum;
			dataProbGivenUnpaired[i] = dataProbGivenUnpaired[i]/sum;
		}
	}
	

	public void transformToAlignment(String gappedseq) {
		int n = gappedseq.length();
		float[] data_a = new float[n];
		float[] dataProbGivenInnerPair_a = new float[n];
		float[] dataProbGivenOuterPair_a = new float[n];
		float[] dataProbGivenUnpaired_a = new float[n];
		
		int cnt = 0; //counts sequence positions
		for(int i = 0; i<n; i++){
			//step alignment positions
			if(MatrixTools.isGap(gappedseq.charAt(i))){
				//if there's a gap, set probabilities to 1
				data_a[i] = -999; //no data for that column
				dataProbGivenInnerPair_a[i] = 1;
				dataProbGivenOuterPair_a[i] = 1;
				dataProbGivenUnpaired_a[i] = 1;
			}
			else{
				data_a[i] = data[cnt];
				dataProbGivenInnerPair_a[i] = dataProbGivenInnerPair[cnt];
				dataProbGivenOuterPair_a[i] = dataProbGivenOuterPair[cnt];
				dataProbGivenUnpaired_a[i] = dataProbGivenUnpaired[cnt];
				cnt++;
			}
			
			//Prevent division by zero by setting 0's to a very small finite number instead
			if(dataProbGivenInnerPair_a[i]==0 && dataProbGivenOuterPair_a[i]==0 && dataProbGivenUnpaired_a[i]==0){
				dataProbGivenInnerPair_a[i]=Float.MIN_VALUE;
				dataProbGivenOuterPair_a[i] =Float.MIN_VALUE;
				dataProbGivenUnpaired_a[i]=Float.MIN_VALUE;
			}
			
			
			//System.out.println("i=" + i + ", data = " + data_a[i] + ": inner pairing="+dataProbGivenInnerPair_a[i] + ", " +
			//		": outer pairing="+dataProbGivenOuterPair_a[i] + ", " +
			//		"unpaired="+dataProbGivenUnpaired_a[i]);
		}
		
		this.data = data_a;
		this.dataProbGivenInnerPair = dataProbGivenInnerPair_a;
		this.dataProbGivenOuterPair = dataProbGivenOuterPair_a;		
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
	}
	
	public void removeColumns(List<Integer> leftoutcolumns){
		Iterator<Integer> iter = leftoutcolumns.iterator();
		int leaveout = 0;
		int from = 0;
		int cnt = 0; //counts position in new thing
		float[] data_a = new float[this.data.length - leftoutcolumns.size()];
		float[] dataProbGivenInnerPair_a = new float[this.data.length - leftoutcolumns.size()];
		float[] dataProbGivenOuterPair_a = new float[this.data.length - leftoutcolumns.size()];
		float[] dataProbGivenUnpaired_a = new float[this.data.length - leftoutcolumns.size()];
		
		while(iter.hasNext()){
			leaveout = iter.next();
			for(int i = from; i<leaveout; i++){
				data_a[cnt] = this.data[i];
				dataProbGivenInnerPair_a[cnt] = dataProbGivenInnerPair[i];
				dataProbGivenOuterPair_a[cnt] = dataProbGivenOuterPair[i];
				dataProbGivenUnpaired_a[cnt] = dataProbGivenUnpaired[i];
				cnt++; 
			}
			from = leaveout+1;
		}
		
		//do the last part part
		for(int i = from; i<data.length; i++){
			data_a[cnt] = this.data[i];
			dataProbGivenInnerPair_a[cnt] = dataProbGivenInnerPair[i];
			dataProbGivenOuterPair_a[cnt] = dataProbGivenOuterPair[i];
			dataProbGivenUnpaired_a[cnt] = dataProbGivenUnpaired[i];
			cnt++; 
		}
		
		
		this.data = data_a;
		this.dataProbGivenInnerPair = dataProbGivenInnerPair_a;
		this.dataProbGivenOuterPair = dataProbGivenOuterPair_a;
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
		//System.out.println("Size of new auxdata data: " + data.length);
	}
	
	
}
