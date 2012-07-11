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


public class SHAPEData implements ExtraData {
	
	private int type = 0; 
	
	//the percent error assumed in each individual SHAPE measurement
	//(this is used to define the interval to calculate probability integrals) 
	
	private float[] data; //contains the actual measurement numbers (data)
	
	//Probabilities must always be given as P(data|unpaired) or P(data|paired).  
	private float[] dataProbGivenPaired; //P(data|model), paired case
	private float[] dataProbGivenUnpaired; //P(data|model), unpaired case
	
	//Unpaired parameter fit
	private final float lamb = 0.681211f;

	private final float k_mixedpair =  0.895341f;
	private final float sigma_mixedpair = 0.0712f;
	private final float mu_mixedpair = 0.00173144f;

	public SHAPEData(){
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
		return dataProbGivenPaired[n];  
	}
	public float getProbabilityGivenOuterPair(int n){
		return dataProbGivenPaired[n];
	}
	public float getProbabilityGivenUnpaired(int n){
		return dataProbGivenUnpaired[n];
	}
	
	private float expPDF(float x, float lamb){
		return (1/lamb) * (float)Math.exp(-x/lamb);
	}
	
	private float gevPDF(float x, float xi, float sigma, float mu){
		float oneoverxi = 1f/xi;
		float oneoversigma = 1f/sigma;
		float xminusmu = x-mu;
		float z = xminusmu*oneoversigma;
		return oneoversigma*(float)Math.exp(-1*Math.pow(1+xi*z,-1*oneoverxi))*
			(float)Math.pow((1+xi*z),-1-oneoverxi); 		
	}
	
	private float calcUnpaired(float f) {
		//Calculates the numbers on the basis of distribution fits
		return expPDF(f,lamb);
	}

	private float calcPaired(float f) {
		// Calculates the numbers on the basis of distribution fits
		return gevPDF(f,k_mixedpair, sigma_mixedpair, mu_mixedpair);
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
		this.dataProbGivenPaired = new float[sequencelength];
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
		//System.out.println("Calculating probs");
		int n = data.length;
		for(int i = 0; i<n; i++){
			//For each data point, calculate the probability of data given paired/unpaired. 
			if(data[i] == -999){
				//ignore 
				dataProbGivenPaired[i] = 1;
				dataProbGivenUnpaired[i] = 1;	
				continue; 
			}

			dataProbGivenPaired[i] = calcPaired(data[i]);
			dataProbGivenUnpaired[i] = calcUnpaired(data[i]);
		}
		
		//For plotting the distributions: 
		//for (int i=0; i<100; i++){
		//	System.out.println(((float)i/100) + " " + calcPaired((float)i/100) + " " + calcUnpaired((float)i/100));
		//}
		//System.exit(0);
	
	}

	public void transformToAlignment(String gappedseq) {
		int n = gappedseq.length();
		float[] data_a = new float[n];
		float[] dataProbGivenPair_a = new float[n];
		float[] dataProbGivenUnpaired_a = new float[n];
		
		int cnt = 0; //counts sequence positions
		for(int i = 0; i<n; i++){
			//step alignment positions
			if(MatrixTools.isGap(gappedseq.charAt(i))){
				//if there's a gap, set probabilities to 1
				data_a[i] = -999; //no data for that column
				dataProbGivenPair_a[i] = 1;
				dataProbGivenUnpaired_a[i] = 1;
			}
			else{
				data_a[i] = data[cnt];
				dataProbGivenPair_a[i] = dataProbGivenPaired[cnt];
				dataProbGivenUnpaired_a[i] = dataProbGivenUnpaired[cnt];
				cnt++;
			}
			
			//Prevent division by zero by setting 0's to a very small finite number instead
			if(dataProbGivenPair_a[i]==0 && dataProbGivenUnpaired_a[i]==0){
				dataProbGivenPair_a[i]=Float.MIN_VALUE;
				dataProbGivenUnpaired_a[i]=Float.MIN_VALUE;
			}
			
			
			//System.out.println("i=" + i + ", data = " + data_a[i] + ": pairing="+dataProbGivenPair_a[i] + 
			//		", unpaired = " + dataProbGivenUnpaired_a[i]);
		}
		
		this.data = data_a;
		this.dataProbGivenPaired = dataProbGivenPair_a;	
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
	}
	
	public void removeColumns(List<Integer> leftoutcolumns){
		Iterator<Integer> iter = leftoutcolumns.iterator();
		int leaveout = 0;
		int from = 0;
		int cnt = 0; //counts position in new thing
		float[] data_a = new float[this.data.length - leftoutcolumns.size()];
		float[] dataProbGivenPair_a = new float[this.data.length - leftoutcolumns.size()];
		float[] dataProbGivenUnpaired_a = new float[this.data.length - leftoutcolumns.size()];
		
		while(iter.hasNext()){
			leaveout = iter.next();
			for(int i = from; i<leaveout; i++){
				data_a[cnt] = this.data[i];
				dataProbGivenPair_a[cnt] = dataProbGivenPaired[i];
				dataProbGivenUnpaired_a[cnt] = dataProbGivenUnpaired[i];
				cnt++; 
			}
			from = leaveout+1;
		}
		
		//do the last part part
		for(int i = from; i<data.length; i++){
			data_a[cnt] = this.data[i];
			dataProbGivenPair_a[cnt] = dataProbGivenPaired[i];
			dataProbGivenUnpaired_a[cnt] = dataProbGivenUnpaired[i];
			cnt++; 
		}
		
		this.data = data_a;
		this.dataProbGivenPaired = dataProbGivenPair_a;
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
		//System.out.println("Size of data after column removal: " + data.length);
	}
	
	private void normalizeProbabilities(){
		System.out.println("Normalizing probs");
		int n = data.length;
		for(int i = 0; i<n; i++){
			//For each data point, normalize the three probabilities 
			//so they are all under 1. 
			//This is so we don't get overflow errors
			//(it's the ratios that matter anyway)
			float sum = dataProbGivenPaired[i] + dataProbGivenUnpaired[i];
			dataProbGivenPaired[i] = dataProbGivenPaired[i]/sum;
			dataProbGivenUnpaired[i] = dataProbGivenUnpaired[i]/sum;
		}
	}
	
	
}
