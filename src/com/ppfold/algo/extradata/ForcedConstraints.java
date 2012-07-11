package com.ppfold.algo.extradata;


import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.ppfold.algo.MatrixTools;

/**
 * Static class to define the model for how to deal with hard constraints.
 * 
 * Input as in Mfold: 
 *
 * force bases i,i+1,...,i+k-1 to be double stranded by entering:
 * F   i   0   k on 1 line in the constraint box.
 *
 * force consecutive base pairs i.j,i+1.j-1, ...,i+k-1.j-k+1 by entering:
 * F   i   j   k on 1 line in the constraint box.
 *
 * force bases i,i+1,...,i+k-1 to be single stranded by entering:
 * P   i   0   k on 1 line in the constraint box.
 *
 * prohibit the consecutive base pairs i.j,i+1.j-1, ...,i+k-1.j-k+1 by entering:
 * P   i   j   k on 1 line in the constraint box.
 *
 * 
 * Types of constraints in data matrix match those in GTfold.
 *
 * data[i][j] for i!=j can have one of the following values: 
 * 0: Nothing is required about a pair i,j in constraints
 * 1: Force pairing between i and j. This implies also that
 *    any non-nested pairings will be prohibited and any 
 *    other pairs i and j would be involved in are also prohibited. 
 * 2: Prohibit pairing between i and j. (Nothing else is done)
 * 
 * data[i][i] can have one of the following values: 
 * 0: Nothing is known about position i at all
 * 3: Position i is forced to be single-stranded. This implies also that
 *    any pairs with i will be prohibited. 
 * 4: Position i is NOT single-stranded, it is forced to be in a pair. 
 * 5: Position i is prohibited from pairing with at least one nucleotide. 
 * 6: Position i is in a forced pair, but the pairing partner is not specified.
 * 
 * Contact distance is implemented as prohibition. 
 * 
 * @author Z.Sukosd
 */


public class ForcedConstraints implements ExtraData {
	
	/*Types of constraints match those in GTfold, for consistent implementation.
	*
	* data[i][j] for i!=j can have one of the following values: 
	* 0: Nothing is required about a pair i,j in constraints
	* 1: Force pairing between i and j. This implies also that
	*    any non-nested pairings will be prohibited and any 
	*    other pairs i and j would be involved in are also prohibited. 
	* 2: Prohibit pairing between i and j. (Nothing else is done)
    * 
	* data[i][i] can have one of the following values: 
	* 0: Nothing is known about position i at all
	* 3: Position i is forced to be single-stranded. This implies also that
	*    any pairs with i will be prohibited. 
	* 4: Position i is NOT single-stranded, it is forced to be in a pair. 
	* 5: Position i is prohibited from pairing with at least one nucleotide. 
	* 6: Position i is in a forced pair, but the pairing partner is not specified.
	*/
	/* 
	*/
	
	int type = 2; 
	
	int [][] data; //This contains the two-dimensional data, but only filled out in one instance. 
	
	/*
	 * These arrays store the raw data: each element is i,j,k, where i,j are the outermost
	 * boundaries and k is the stem length.
	 */
	ArrayList<int[]> forceArray = new ArrayList<int[]>();
	ArrayList<int[]> prohibitArray = new ArrayList<int[]>();
	int contactDistance = -1; //Maximum contact distance; -1 if not set 
	
	public ForcedConstraints(){}

	public int getType(){
		return type; 
	}

	
	public ForcedConstraints combinedForcedConstraints(int n, ArrayList<ForcedConstraints> constraints){
		//This returns a specific instance of forced constraints with an empty matrix.
		//This is important because not every ForcedConstraints will have the data matrix
		//initialized, as it fills too much in memory. Only one (the "main") ForcedConstraints
		//will have the datamatrix initialized. 
		ForcedConstraints fcons = new ForcedConstraints();
		fcons.data = new int[n][n]; //Set all data to 0's
		for(int i = 0; i<n; i++){
			for(int j = 0; j<n; j++){
				fcons.data[i][j] =0;
 			}
		}

		for(ForcedConstraints constraint:constraints){
			//System.out.println(toString(forceArray));
			//System.out.println(toString(prohibitArray));
			
			fcons.forceArray.addAll(constraint.forceArray);
			fcons.prohibitArray.addAll(constraint.prohibitArray);

			//replace with minimum contact distance  (to begin with it's -1)
			if(fcons.contactDistance==-1 || (fcons.contactDistance > -1 && constraint.contactDistance<fcons.contactDistance)){
				fcons.contactDistance = constraint.contactDistance;
			}
		}
		
		if(fcons.contactDistance > -1){
			System.out.println("Prevailing contact distance is: " + fcons.contactDistance);
		}

		//Now do the conversion to the table
		//Handle prohibiting and single-stranded 
		for(int[] cons:fcons.prohibitArray){
			//Smallest coordinate should come first 
			if(cons[1] > -1 && cons[0] > cons[1]){
				int tmp = cons[0];
				cons[0] = cons[1];
				cons[1] = tmp; 
			}
			for(int k = 1; k <= cons[2]; k++){
				int itarget = cons[0]+k-1;
				int jtarget = cons[1]-(k-1);

				if(cons[1] == -1){
					//force single-stranded 
					//(No check needed as 3 is stronger than 5)
					fcons.data[itarget][itarget] = 3;
					//prohibit all pairs with that base
					for(int i = 0; i<n; i++){
						fcons.data[i][itarget] = 2;
						fcons.data[itarget][i] = 2;
					}
				}
				else{
					//prohibit the pair
					fcons.data[itarget][jtarget] = 2;
					fcons.data[jtarget][itarget] = 2;
					//If it's not forced to pair or stay unpaired, mark that these two nucleotides are involved in a prohibited pair
					//Need to check as 5 is weaker than 4 or 3
					if(fcons.data[itarget][jtarget]!=4&&fcons.data[itarget][jtarget]!=3){
						fcons.data[itarget][itarget] = 5;
						fcons.data[jtarget][jtarget] = 5;
					}
				}
			}
		}
		
		//prohibit pairs too close
		if(fcons.contactDistance > -1){
			for (int i = 0; i < n; i++) {
				for (int j = i+1; j < n; j++) {
					if (j - i > fcons.contactDistance) {
						fcons.data[i][j] = 2;
						fcons.data[j][i] = 2;
						//System.out.println(i + " and " + (i+j) + " prohibited (too close)");
					}
					if(fcons.data[i][i]!=4&&fcons.data[i][i]!=3){
						//If it's not forced to pair or stay unpaired, mark that it's prohibited from pairing.
						//Need to check as 5 is weaker than 4 or 3
						fcons.data[i][i] = 5;
					}
					if(fcons.data[j][j]!=4&&fcons.data[j][j]!=3){
						//If it's not forced to pair or stay unpaired, mark that it's prohibited from pairing.
						//Need to check as 5 is weaker than 4 or 3
						fcons.data[j][j] = 5;
					}
				}
			}
		}

		//Handle forcing
		for(int[] cons:fcons.forceArray){
			//Smallest coordinate should come first
			if(cons[1] > -1 && cons[0] > cons[1]){
				int tmp = cons[0];
				cons[0] = cons[1];
				cons[1] = tmp; 
			}			
			for(int k = 1; k<=cons[2]; k++){
				//force pairing
				int itarget = cons[0]+k-1;
				int jtarget = cons[1]-(k-1);

				if(cons[1] != -1){
					//prohibit all pairs not-nested with respect to this one 
					//(including the ones which include either of the bases)
					for(int i = 0; i<itarget; i++){
						for(int j = itarget; j<=jtarget; j++){
							fcons.data[i][j]=2;
							fcons.data[j][i]=2;
						}
					}
					for(int i = itarget+1; i<jtarget; i++){
						for(int j = jtarget+1; j<n; j++){
							fcons.data[i][j]=2;
							fcons.data[j][i]=2;
						}
					}
					//prohibit all remaining pairs with the pairing bases 
					//(inside enclosed region)
					for(int i = itarget+1; i<jtarget; i++){
						fcons.data[itarget][i] = 2;
						fcons.data[i][itarget] = 2;
						fcons.data[jtarget][i] = 2;
						fcons.data[i][jtarget] = 2;
					}

					//mark that these two nucleotides are involved in a constrained pair
					//to avoid searching in O(N^2) time
					//(No check needed as 4 is stronger than 5)
					fcons.data[itarget][itarget] = 4;
					fcons.data[jtarget][jtarget] = 4;
					//force pairing 
					fcons.data[itarget][jtarget] = 1;
					fcons.data[jtarget][itarget] = 1;
				}
				else{
					//force pairing, partner is not known 
					//Check that it's not forced to pair (it shouldn't be forced to stay unpaired simultaneously)
					if(fcons.data[itarget][itarget]!=4){
						fcons.data[itarget][itarget] = 6;
					}
					//prohibit all pairs with that base
				}
			}
		}
		//MatrixTools.print(fcons.data);
		return fcons;
	}
	public boolean isEmpty(int i)
	{
		//System.out.println(toString(forceArray));
		//System.out.println(toString(prohibitArray));
		for(int[] cons:prohibitArray){
			if(cons[0]==i || cons[1]==i){
				//System.out.println("data for col " + i + " not empty ");
				return false;
			}
		}
		for(int[] cons:forceArray){
			if(cons[0]==i || cons[1]==i){
				//System.out.println("data for col " + i + " not empty ");
				return false;
			}
		}
		//System.out.println("data for col " + i + " is empty ");
		return true;
	}

	public float getProbabilityGivenOuterPaired(int position1, int position2) {
		return canPair(position1, position2);
	}
	public float getProbabilityGivenInnerPaired(int position1, int position2) {
		return canPair(position1, position2);
	}

	public float getProbabilityGivenUnpaired(int n){
		return canSs(n);
	}

	private float canPair(int i, int j){
		//can't pair if i,j is prohibited, or i or j are forced to be single-stranded
		if(data[i][j]==2||data[i][i]==3||data[j][j]==3){return Float.MIN_VALUE;}
		//can't pair if i or j are forced to pair with something other than each other
		if((data[i][i]==4||data[j][j]==4)&&data[i][j]!=1){return Float.MIN_VALUE;}
		return 1f; 
	}

	private float canSs(int i){
		if(data[i][i]!=4&&data[i][i]!=6){return 1f;}
		else{return Float.MIN_VALUE;}
	}

	/**
	 * Reads SHAPE distribution data from the file specified 
	 * @throws IOException 
	 */

	public void importData(String filename, int sequencelength) throws Exception{
		//read the data 
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			System.err.println("Constraint sequence input file " + filename + " could not be read!");
			throw new IOException(e);
		}
		readData_toStream(stream, sequencelength);	
	}
	
	public void readData_toStream(BufferedInputStream stream, int sequencelength) throws Exception{
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
				boolean [] readdata_letter =  new boolean[data_size];
				// Create a pattern to match different kinds of separators
		        Pattern p = Pattern.compile("[,\\s]+");
		        //DO NOT ignore first line... 
				for(int i = 0; i<lines.length; i++){
					line = lines[i];
					String splitline[] = p.split(line.trim());
					if(splitline.length == 4){
						//True if prohibiting, False if forcing
						readdata_letter[i] = (Character.toLowerCase(splitline[0].charAt(0))=='p')?true:false; 
						int [] readdata = new int[3];
						readdata[0] = Integer.valueOf(splitline[1])-1; //start of pairing
						//note data files are numbered from 1; Java numbers from 0
						readdata[1] = Integer.valueOf(splitline[2])-1; //end of pairing
						readdata[2] = Integer.valueOf(splitline[3]); //length of pairing
						if(!readdata_letter[i]){
							forceArray.add(readdata);
						}
						else{
							prohibitArray.add(readdata);
						}
					}
				}
			}
			catch(Exception e){
				System.err.println("An exception occured while attempting to read or interpret the constraint data. ");
				throw new Exception(e);
			}
		}
		else{
			System.err.println("BufferedStream was null, ignoring...");
			//The stream is null if ONLY contact distance is given. (or maybe not even that)
			//In this case, do nothing, the object will just be empty.
		}
	}
	
	public void setContactDistance(int n){
		contactDistance = n;
	}
	public int getContactDistance(){
		return contactDistance;
	}
	
	public void transformToAlignment(String gappedseq) {
		//System.out.println("before transformtoalignment");
		//System.out.println(toString(forceArray));
		//System.out.println(toString(prohibitArray));
		if(gappedseq == null){
			//we just have a contact distance, no array data.
			return;
		}
		int n = gappedseq.length();
		int[] nrofgaps = new int[n];
		int gapcounter = 0; //counts gaps
		for(int i = 0; i<n; i++){
			//step alignment positions
			if(MatrixTools.isGap(gappedseq.charAt(i))){
				gapcounter++;
			}
			nrofgaps[i] = gapcounter;
		}		
		//Transform with the gaps
		for(int[] cons:prohibitArray){
			cons[0] = cons[0] + nrofgaps[cons[0]];
			if(cons[1]>-1){ //this is to make sure P i 0 k doesn't get altered
				cons[1] = cons[1] + nrofgaps[cons[1]];
			}
		}
		for(int[] cons:forceArray){ //this is to make sure F i 0 k doesn't get altered
			cons[0] = cons[0] + nrofgaps[cons[0]];
			if(cons[1]>-1){
				cons[1] = cons[1] + nrofgaps[cons[1]];
			}
		}
		//System.out.println("after transformtoalignment");
		//System.out.println(toString(forceArray));
		//System.out.println(toString(prohibitArray));
	}
	
	public void removeColumns(List<Integer> leftoutcolumns) throws Exception{
		//System.out.println("before removing columns");
		//System.out.println(toString(forceArray));
		//System.out.println(toString(prohibitArray));
		if(forceArray.size()==0&&prohibitArray.size()==0){
			//In this case we don't actually have any data except maybe contact distance, so do nothing 
			return;
		}
		@SuppressWarnings("unchecked")
		int maxpos = findMaxPos(forceArray,prohibitArray);
		int [] remColCnt = new int[maxpos+1];
		Iterator<Integer> iter = leftoutcolumns.iterator();
		int current = 0;
		int cnt = 0; 
		while(iter.hasNext()){
			int col = iter.next();
			for(int i = current; i<=col; i++){
				if(i == col){
					cnt ++; //one more column removed
					remColCnt[i] = -1; //Remove this column 
					current = i+1; 
					break;
				}
				else{
					remColCnt[i] = cnt; 
				}
			}
		}
		//fill last bit
		for(int i = current; i<=maxpos; i++){
			remColCnt[i] = cnt; 
		}
		
		//Transform with the removed columns
		//In theory, there can't be any columns that are kept that have associated constraints. 
		//If so, throw an exception... 
		for(int[] cons:prohibitArray){
			if(remColCnt[cons[0]]!=-1){
				cons[0] = cons[0] - remColCnt[cons[0]];
			}
			else{
				throw new Exception("Column with data was trying to be removed!");
			}
			if(cons[1]>-1){
				if(remColCnt[cons[1]]!=-1){
					cons[1] = cons[1] - remColCnt[cons[1]];
				}
				else{
					throw new Exception("Column with data was trying to be removed!");
				}
			}
		}
		for(int[] cons:forceArray){
			if(remColCnt[cons[0]]!=-1){
				cons[0] = cons[0] - remColCnt[cons[0]];
			}
			else{
				throw new Exception("Column with data was trying to be removed!");
			}
			if(cons[1]>-1){
				if(remColCnt[cons[1]]!=-1){
					cons[1] = cons[1] - remColCnt[cons[1]];
				}
				else{
					throw new Exception("Column with data was trying to be removed!");
				}
			}
		}
		//System.out.println("after removing columns");
		//System.out.println(toString(forceArray));
		//System.out.println(toString(prohibitArray));
		
		
	}

	private static String toString(ArrayList<int[]> inputlist){
		String result = "";
		for(int[] input:inputlist ){
		result = result + "[ ";
		for(int i:input){
			result = result +i + " ";
		}
		result = result + "] \n";
		}
		return result;
	}
	
	private int findMaxPos(ArrayList<int[]> ... arrays){
		int res = 0;
		for(ArrayList<int[]> array : arrays){
			for(int[] cons : array){
				if(res < cons[0]){
					res = cons[0];
				}
				if(res < cons[1]){
					res = cons[1];
				}
			}
		}
		return res;
	}
	
	
	
}
