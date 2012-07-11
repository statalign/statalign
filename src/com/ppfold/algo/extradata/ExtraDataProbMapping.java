package com.ppfold.algo.extradata;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.ppfold.algo.MatrixTools;

public class ExtraDataProbMapping implements ExtraData {
	//This is a special type of ExtraData, where the numbers in the input file are actually
	//the probability densities of some observation given paired/unpaired. It is assumed the
	//user has calculated these numbers on their own according to some model. 
	
	private int type = 1; 
	
	//Probabilities must always be given as P(data|unpaired) or P(data|paired).  
	private float[] dataProbGivenPaired; //P(data|model), unpaired case 
	private float[] dataProbGivenUnpaired; //P(data|model), paired case
	
	public int getType(){
		return type; 
	}
	public float getProbabilityGivenOuterPaired(int position1, int position2) {
		return dataProbGivenPaired[position1]*dataProbGivenPaired[position2];
	}
	public float getProbabilityGivenInnerPaired(int position1, int position2) {
		return dataProbGivenPaired[position1]*dataProbGivenPaired[position2];
	}

	public float getProbabilityGivenUnpaired(int position) {
		return dataProbGivenUnpaired[position];
	}


	public void importData(String filename, int sequencelength) throws Exception{
		//read the SHAPE data 
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			System.err.println("Extra data input file " + filename + " could not be read!");
			throw new IOException(e);
		}
		readData_toStream(stream, sequencelength);	
	}
	

	public void readData_toStream(BufferedInputStream stream, int sequencelength) throws Exception{
		//create the data containers
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
				float [] readdata_data1 = new float[data_size];
				float[] readdata_data2 = new float[data_size];
				// Create a pattern to match different kinds of separators
		        Pattern p = Pattern.compile("[,\\s]+");
		        //DO NOT ignore first line... 
				for(int i = 0; i<lines.length; i++){
					line = lines[i];
					String splitline[] = p.split(line.trim());
					if(splitline.length == 3){
						readdata_index[i] = Integer.valueOf(splitline[0])-1; //data files are numbered from 1; Java numbers from 0
						readdata_data1[i] = Float.valueOf(splitline[1]);
						readdata_data2[i] = Float.valueOf(splitline[2]);
					}
				}

				//initialize all data values to -999 and all probability values to 1
				for(int i = 0; i<sequencelength; i++){
					dataProbGivenPaired[i] = 1;
					dataProbGivenUnpaired[i] = 1;
				}
				//set the appropriate data values to the read ones
				for(int i = 0; i<data_size; i++){
					//
					dataProbGivenPaired[readdata_index[i]] = readdata_data1[i];
					dataProbGivenUnpaired[readdata_index[i]] = readdata_data2[i];
				}
				
			}
			catch(Exception e){
				System.err.println("An exception occured while attempting to read or interpret the data. ");
				throw new Exception(e);
			}
		}
		else{
			System.err.println("Input stream was null, the data could not be loaded.");
		}
	}

	public void transformToAlignment(String gappedseq) {
		int n = gappedseq.length();
		float[] dataProbGivenPaired_a = new float[n];
		float[] dataProbGivenUnpaired_a = new float[n];
		
		int cnt = 0; //counts sequence positions
		for(int i = 0; i<n; i++){
			//step alignment positions
			if(MatrixTools.isGap(gappedseq.charAt(i))){
				//if there's a gap, set probabilities to 1
				dataProbGivenPaired_a[i] = 1;
				dataProbGivenUnpaired_a[i] = 1;
			}
			else{
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
			
			
			//System.out.println(i + ": pairing="+dataProbGivenPaired_a[i] + ", " +
			//		"unpaired="+dataProbGivenUnpaired_a[i]);
		}
		
		this.dataProbGivenPaired = dataProbGivenPaired_a;
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
		
	}

	public void removeColumns(List<Integer> leftoutcolumns){
		Iterator<Integer> iter = leftoutcolumns.iterator();
		int leaveout = 0;
		int from = 0;
		int cnt = 0; //counts position in new thing
		int length = this.dataProbGivenPaired.length; 
		float[] dataProbGivenPaired_a = new float[length - leftoutcolumns.size()];
		float[] dataProbGivenUnpaired_a = new float[length - leftoutcolumns.size()];
		
		while(iter.hasNext()){
			leaveout = iter.next();
			for(int i = from; i<leaveout; i++){
				dataProbGivenPaired_a[cnt] = dataProbGivenPaired[i];
				dataProbGivenUnpaired_a[cnt] = dataProbGivenUnpaired[i];
				cnt++; 
			}
			from = leaveout+1;
		}
		
		//do the last part part
		for(int i = from; i<length; i++){
			dataProbGivenPaired_a[cnt] = dataProbGivenPaired[i];
			dataProbGivenUnpaired_a[cnt] = dataProbGivenUnpaired[i];
			cnt++; 
		}
		
		this.dataProbGivenPaired = dataProbGivenPaired_a;
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
		//System.out.println("Size of new aux data: " + dataProbGivenPaired.length);
	}
	
	
	public boolean isEmpty(int i)
	{
		return (dataProbGivenPaired[i]==1&&dataProbGivenUnpaired[i]==1);
	}

}
