package com.ppfold.algo.extradata;


import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.ppfold.algo.MatrixTools;



/**
 * Static class to define the model for how to deal with SHAPE data as well as container for the actual data 
 * 
 * @author Z.Sukosd
 */


public class ExtraDataBars implements ExtraData {
	
	//the percent error assumed in each individual SHAPE measurement
	//(this is used to define the interval to calculate probability integrals) 
	
	private int type = 0; 
	
	private float[] data; //contains the actual measurement numbers (data)
	
	private float[] distributionPaired; //contains the probability density function, paired case (model)
	
	private float[] distributionUnpaired; //contains the probability density function, unpaired case (model) 
	private float[] distributionLimits; //contains the LOWER limits for the probability densities (model)
	
	//Probabilities must always be given as P(data|unpaired) or P(data|paired).  
	private float[] dataProbGivenPaired; //P(data|model), unpaired case 
	private float[] dataProbGivenUnpaired; //P(data|model), paired case

	public ExtraDataBars(){}
	
	public int getType(){
		return type; 
	}

	
	public boolean isEmpty(int i)
	{
		return data[i]==-999;
	}
	
	private void setDistribution(float[] unpaired, float[] paired, float[] limits){
		distributionPaired = paired;
		distributionUnpaired = unpaired;
		distributionLimits = limits;
	}

	public float getProbabilityGivenOuterPaired(int position1, int position2) {
		return getProbabilityGivenPaired(position1)*getProbabilityGivenPaired(position2);
	}
	public float getProbabilityGivenInnerPaired(int position1, int position2) {
		return getProbabilityGivenPaired(position1)*getProbabilityGivenPaired(position2);
	}
	
	public float[] getDistributionPaired() {
		return distributionPaired;
	}
	public float[] getDistributionUnpaired() {
		return distributionUnpaired;
	}
	public float[] getDistributionLimits() {
		return distributionLimits;
	}
	

	public float getProbabilityGivenPaired(int n){
		return dataProbGivenPaired[n];  
	}
	
	public float getProbabilityGivenInnerPair(int n){
		return getProbabilityGivenPaired(n);  
	}
	
	public float getProbabilityGivenOuterPair(int n){
		return getProbabilityGivenPaired(n);  
	}
	
	public float getProbabilityGivenUnpaired(int n){
		return dataProbGivenUnpaired[n];
	}
	
	/**
	 * Reads SHAPE distribution data from the file specified 
	 * @throws IOException 
	 */

	public static ExtraDataBars readDistTable(String filename) throws IOException{
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			System.err.println("SHAPE input file " + filename + " could not be read!");
			throw new IOException(e);
		}
		return readDistTable_toStream(stream);
	}
	
	public static ExtraDataBars readDistTable_toStream(BufferedInputStream stream){
		//the CLC import plugin uses the following line:  
		//public ClcObject[] doImport(BufferedInputStream stream, String name, Activity activity) throws IOException, ParseException, PersistenceException {
		//the same stream will be passed here to ease the import process later. 
		ExtraDataBars result = null;
		if(stream!=null){
			try{
				int l = stream.available();
				byte[] bytes = new byte[l];
				stream.read(bytes);
				stream.close();
				String data_string = new String(bytes);
				String[] lines = data_string.split("\n");
				//First line is text 
				int data_size = lines.length-1;
				//Create SHAPEdata object 
				result = new ExtraDataBars();
				float [] paired = new float[data_size];
				float [] unpaired = new float[data_size];
				float [] limits = new float[data_size];
				// Create a pattern to match different kinds of separators
		        Pattern p = Pattern.compile("[,\\s]+");
		        //Ignore first line... 
				for(int i = 1; i<lines.length; i++){
					String line = lines[i];
					String splitline[] = p.split(line.trim());
					limits[i-1] = Float.valueOf(splitline[0]);
					paired[i-1] = Float.valueOf(splitline[1]);
					unpaired[i-1] = Float.valueOf(splitline[2]);
				}
				result.setDistribution(unpaired,paired,limits);
			}
			catch(Exception e){
				System.err.println("An exception occured while attempting to read the SHAPE data.");
			}

		}
		else{
			System.err.println("Input stream was null, SHAPE data could not be loaded.");
		}
		return result;
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
		        boolean validdata = false;
		        //DO NOT ignore first line...
				for(int i = 0; i<lines.length; i++){
					line = lines[i];
					String splitline[] = p.split(line.trim());
					if(splitline.length == 2){
						validdata = true;
						readdata_index[i] = Integer.valueOf(splitline[0])-1; //SHAPE data files are numbered from 1; Java numbers from 0
						readdata_data[i] = Float.valueOf(splitline[1]);
					}
				}

				//initialize all data values to -999
				for(int i = 0; i<this.data.length; i++){
					data[i] = -999;
				}
				//set the appropriate data values to the read ones
				for(int i = 0; validdata&&i<data_size; i++){
					data[readdata_index[i]] = readdata_data[i];
				}
				
				//calculate probabilities 
				try{
					calcProbabilities();
				}
				catch(Exception e){
					System.out.println("Probabilities for the data could not be calculated.");
					throw new Exception(e);
				}
			}
			catch(Exception e){
				System.err.println("An exception occured while attempting to read or interpret the SHAPE sequence data. ");
				e.printStackTrace();
				throw new Exception(e);
			}
		}
		else{
			System.err.println("Input stream was null, SHAPE sequence data could not be loaded.");
		}
	}
	
	private void calcProbabilities(){
		int n = data.length;
		int distlength = distributionLimits.length;
		for(int i = 0; i<n; i++){
			//For each data point, calculate the probability of data given paired/unpaired. 
			if(data[i] == -999){
				//ignore 
				dataProbGivenPaired[i] = 1;
				dataProbGivenUnpaired[i] = 1;
				continue; 
			}
			int index = 0;
			while(index!=distlength && !(data[i] < distributionLimits[index])){
				index++;
			}
			index--;

			dataProbGivenPaired[i] = distributionPaired[index];
			dataProbGivenUnpaired[i] = distributionUnpaired[index];				

			//System.out.println(data[i] + ": pairing="+dataProbGivenPaired[i] + ", " +
			//		"unpaired="+dataProbGivenUnpaired[i]);
		}
	}

	
	public void transformToAlignment(String gappedseq) {
		int n = gappedseq.length();
		float[] data_a = new float[n];
		float[] dataProbGivenPaired_a = new float[n];
		float[] dataProbGivenUnpaired_a = new float[n];
		
		int cnt = 0; //counts sequence positions
		for(int i = 0; i<n; i++){
			//step alignment positions
			if(MatrixTools.isGap(gappedseq.charAt(i))){
				//if there's a gap, set probabilities to 1
				data_a[i] = -999; //no data for that column
				dataProbGivenPaired_a[i] = 1;
				dataProbGivenUnpaired_a[i] = 1;
			}
			else{
				data_a[i] = data[cnt];
				dataProbGivenPaired_a[i] = dataProbGivenPaired[cnt];
				dataProbGivenUnpaired_a[i] = dataProbGivenUnpaired[cnt];
				cnt++;
			}
			
			//Prevent weird results by setting 0's to a very small finite number instead
			if(dataProbGivenPaired_a[i]==0){
				dataProbGivenPaired_a[i]=Float.MIN_VALUE;
			}
			if(dataProbGivenUnpaired_a[i]==0){
				dataProbGivenUnpaired_a[i]=Float.MIN_VALUE;
			}
			
			//System.out.println(data_a[i] + ": pairing="+dataProbGivenPaired_a[i] + ", " +
			//		"unpaired="+dataProbGivenUnpaired_a[i]);
		}
		
		this.data = data_a;
		this.dataProbGivenPaired = dataProbGivenPaired_a;
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
	}
	
	public void removeColumns(List<Integer> leftoutcolumns){
		Iterator<Integer> iter = leftoutcolumns.iterator();
		int leaveout = 0;
		int from = 0;
		int cnt = 0; //counts position in new thing
		float[] data_a = new float[this.data.length - leftoutcolumns.size()];
		float[] dataProbGivenPaired_a = new float[this.data.length - leftoutcolumns.size()];
		float[] dataProbGivenUnpaired_a = new float[this.data.length - leftoutcolumns.size()];
		
		while(iter.hasNext()){
			leaveout = iter.next();
			for(int i = from; i<leaveout; i++){
				data_a[cnt] = this.data[i];
				dataProbGivenPaired_a[cnt] = dataProbGivenPaired[i];
				dataProbGivenUnpaired_a[cnt] = dataProbGivenUnpaired[i];
				cnt++; 
			}
			from = leaveout+1;
		}
		
		//do the last part part
		for(int i = from; i<data.length; i++){
			data_a[cnt] = this.data[i];
			dataProbGivenPaired_a[cnt] = dataProbGivenPaired[i];
			dataProbGivenUnpaired_a[cnt] = dataProbGivenUnpaired[i];
			cnt++; 
		}
		
		
		this.data = data_a;
		this.dataProbGivenPaired = dataProbGivenPaired_a;
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
		//System.out.println("Size of new auxdata data: " + data.length);
	}

	
	
}
