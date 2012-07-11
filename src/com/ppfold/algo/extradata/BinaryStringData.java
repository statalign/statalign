package com.ppfold.algo.extradata;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.ppfold.algo.MatrixTools;

public class BinaryStringData implements ExtraData {

	private int type = 0; 
	private int[] data; //contains the actual measurement numbers (data)
	//Probabilities must always be given as P(data|unpaired) or P(data|paired).  
	private float[] dataProbGivenPaired; //P(data|model), unpaired case 
	private float[] dataProbGivenUnpaired; //P(data|model), paired case
	
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
	
	public float getProbabilityGivenInnerPair(int position) {
		return dataProbGivenPaired[position];
	}

	public float getProbabilityGivenOuterPair(int position) {
		return dataProbGivenPaired[position];
	}

	public float getProbabilityGivenUnpaired(int position) {
		return dataProbGivenUnpaired[position];
	}

	public void importData(String filename, int sequencelength) throws Exception{
		//read the data 
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			System.err.println("Binary string input file " + filename + " could not be read!");
			throw new Exception(e);
		}
		readData_toStream(stream, sequencelength);	
	}
	
	public void readData_toStream(BufferedInputStream stream, int sequencelength) throws Exception{
		//the CLC import plugin uses the following line:  
		//public ClcObject[] doImport(BufferedInputStream stream, String name, Activity activity) throws IOException, ParseException, PersistenceException {
		//the same stream will be passed here to ease the import process later.
		
		//create the data containers
		this.data = new int[sequencelength];
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
				int [] readdata_data = new int[data_size];
				// Create a pattern to match different kinds of separators
		        Pattern p = Pattern.compile("[,\\s]+");
		        //DO NOT ignore first line... 
				for(int i = 0; i<lines.length; i++){
					line = lines[i];
					String splitline[] = p.split(line.trim());
					if(splitline.length == 2){
						readdata_index[i] = Integer.valueOf(splitline[0])-1; //SHAPE data files are numbered from 1; Java numbers from 0
						readdata_data[i] = Integer.valueOf(splitline[1]);
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
			}
			catch(Exception e){
				System.err.println("An exception occured while attempting to read or interpret the binary string sequence data. ");
				throw new Exception(e);
			}
		}
		else{
			System.err.println("Input stream was null, binary string sequence data could not be loaded.");
		}
	}
	
	private void calcProbabilities(){
		int n = data.length;
		for(int i = 0; i<n; i++){
			//For each data point, calculate the probability of data given paired/unpaired. 
			//If no data, set both to 1 (so other info fills them in)
			if(data[i] == -999){
				//ignore 
				dataProbGivenPaired[i] = 1;
				dataProbGivenUnpaired[i] = 1;
				continue; 
			}

			dataProbGivenPaired[i] = data[i]==1?1:Float.MIN_VALUE;
			dataProbGivenUnpaired[i] = data[i]==1?Float.MIN_VALUE:1;				
		}
	}
	
	public void transformToAlignment(String gappedseq) {
		int n = gappedseq.length();
		int[] data_a = new int[n];
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
			
			//Prevent division by zero by setting 0's to a very small finite number instead
			if(dataProbGivenPaired_a[i]==0 && dataProbGivenUnpaired_a[i]==0){
				dataProbGivenPaired_a[i]=Float.MIN_VALUE;
				dataProbGivenUnpaired_a[i]=Float.MIN_VALUE;
			}
			
			
			//System.out.println((i+1) + ": data: " + data_a[i] + ", values: pairing="+dataProbGivenPaired_a[i] + ", " +
			//		"unpaired="+dataProbGivenUnpaired_a[i]);
		}
		
		this.data = data_a;
		this.dataProbGivenPaired = dataProbGivenPaired_a;
		this.dataProbGivenUnpaired = dataProbGivenUnpaired_a;	
	}
	
	public void removeColumns(List<Integer> leftoutcolumns) {
		Iterator<Integer> iter = leftoutcolumns.iterator();
		int leaveout = 0;
		int from = 0;
		int cnt = 0; //counts position in new thing
		int[] data_a = new int[this.data.length - leftoutcolumns.size()];
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
	
}
	
	
}
