package com.ppfold.main;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;

import javax.swing.JDialog;

import com.ppfold.algo.MatrixTools;
import com.ppfold.algo.Node;
import com.ppfold.algo.Tree;
import com.ppfold.algo.extradata.ExtraDataBars;

/**
 * Executes SOME checks of the input data. It is not guaranteed to find all possible errors.
 * 
 * @author Z.Sukosd
 */

public class CheckAllData {
	static boolean alignmentOK;
	static boolean distOK; 
	static boolean dataOK;
	static boolean fileTooLong;
	static boolean nothingfailed;
	static String reason;
	public static final String SUCCESS_TEXT = "THERE WERE FAILED TESTS! Details below.\n\n";
	public static final String FAILURE_TEXT = "EVERYTHING SEEMS TO BE OK! Details below.\n\n";
	private static final String gapstring = "[-.]";
	
	static String checkData() throws Exception {
		alignmentOK = true;
		distOK = true;
		dataOK = true; 
		nothingfailed = true;
		String testresults = "";
		reason = "";
		//Check alignment 
		testresults += "Alignment check result: \n" ;
		fileTooLong=false;
		Alignment ali = attemptToLoadAlignment();
		
		boolean alignmentCanBeOpened = ali==null?false:true;
		testresults += "- File exists and can be opened... ";
		if(fileTooLong){
			testresults += "FAILED (cause: file too large)\n";
		}
		else{
			testresults += alignmentCanBeOpened?"OK\n":"FAILED\n";
		}
	
		if(alignmentCanBeOpened){
			testresults += "- Sequences match names... ";
			testresults += namesMatchSequences(ali)?"OK\n":"FAILED\n";
			testresults += "- No illegal symbols... ";
			String symbol = noIllegalSymbols(ali);
			testresults += symbol==null?"OK\n":"FAILED (cause: " + symbol + ")\n";
			testresults += "- Gapped sequences have the same length... ";
			int nr = lengthOK(ali);
			if(nr>0){testresults += "FAILED (cause: sequence nr. " + (nr+1) + ")\n";}
			else if(nr==-1){testresults += "FAILED (cause: general format error)\n";}
			else{testresults += "OK\n";}
			
			testresults += "- Memory requirements met... ";
			long required = memRequired(ali);
			long available = Runtime.getRuntime().maxMemory()/(1024*1024); 
			if(required<available){
				testresults += "OK (required estimate: " + required + " MB, available: " + Runtime.getRuntime().maxMemory()/(1024*1024) + " MB)\n";
			}
			else{
				nothingfailed = false;
				testresults += "FAILED (required estimate: " + required + " MB, available: " + Runtime.getRuntime().maxMemory()/(1024*1024) + " MB)\n"; 
			}
				
		}
	
		testresults += "\n";
		testresults += "Tree check result: \n";
		
		boolean treeCanBeOpened = false;
		Tree tree = null;
		if(PPfoldMain.treefilename==null){
			testresults += "- Tree not provided \n";
		}
		else{
			tree = attemptToLoadTree();
			treeCanBeOpened = tree==null?false:true;
			testresults += "- File exists and can be parsed... ";
			if(fileTooLong){
				testresults += "FAILED (cause: file too large)\n";
			}
			else{
				testresults += treeCanBeOpened?"OK\n":"FAILED\n";
			}		}

		if(treeCanBeOpened&&alignmentOK){
			testresults += "- Tree names match alignment names... ";
			String results = namesMatchTreeAlignment(ali,tree);
			testresults += results==null?"OK\n":"FAILED (cause: node named " + results +")\n";
		}

		testresults += "\n";
		testresults += "Output folder check result: \n";
		boolean outputExists = attemptToAccessOutput(); 
		testresults += "- Folder exists... ";
		testresults += outputExists?"OK\n":"FAILED\n";

		if(outputExists){
			boolean outputFolderWriteable = attemptToWriteOutput();
			testresults += "- Output folder writeable... ";
			testresults += outputFolderWriteable?"OK\n":"FAILED\n";
		}

		testresults += "\n";
		testresults += "Data check result:\n";
		if(PPfoldMain.datainfo.size()==0){
			testresults += "- Data not provided. \n";
		}
		else{
			for(DataInfo data:PPfoldMain.datainfo){
				dataOK = true;
				if(data.getType()==0){
					boolean dataexists;
					ExtraDataBars dist;
					boolean distheader;
					boolean distprobs;
					int seqID;
					distOK = true;
					testresults += "- Data identifier: " + data.getiD() + "\n";

					testresults += "  - Distribution file can be opened: ";

					dist = attemptToLoadDataDist(data);
					testresults += dist!=null?"OK\n":"FAILED\n";

					if(dist!=null){
						testresults += "  - Distribution file has the correct header: ";
						distheader = checkDistHeader(data);
						testresults += distheader?"OK\n":"FAILED\n";

						testresults += "  - Distribution probabilities add to 1 and limits are increasing: "; 
						distprobs = checkDistProbs(dist);
						testresults += distprobs?"OK\n":("FAILED (cause: "+ reason + ")\n");	
					}

					dataexists = attemptToLoadData(data);
					testresults += "  - Data file can be opened: "; 
					testresults += dataexists?"OK\n":"FAILED\n";

					int[] indices = null; 
					if(dataexists){
						indices = checkBarDataFormat(data);
						testresults += "  - Data has right format: "; 
						testresults += indices!=null?"OK\n":"FAILED\n";					
					}

					if(alignmentOK){
						seqID = sequenceFound(ali, data);
						testresults += "  - Sequence name found in alignment: "; 
						testresults += seqID!=-1?"OK\n":"FAILED\n";

						if(seqID!=-1&&dataOK&&dataexists){
							testresults += "  - Data not longer than sequence: "; 
							String newseq = new String(ali.getSequences().get(seqID));
							newseq.replaceAll(gapstring, "");
							testresults += sequenceMatchesBarData(newseq.length(), indices)?"OK\n":"FAILED\n";	
						}

					}
				}
				else if(data.getType()==1){
					int seqID; 
					boolean dataexists = false; 
					testresults += "- Data identifier: " + data.getiD() + "\n";
					
					dataexists = attemptToLoadData(data);
					testresults += "  - Data file can be opened: "; 
					testresults += dataexists?"OK\n":"FAILED\n";

					int[] indices = null; 
					if(dataexists){
						indices = checkProbDataFormat(data);
						testresults += "  - Data has right format: "; 
						testresults += indices!=null?"OK\n":"FAILED\n";					
					}

					if(alignmentOK){
						seqID = sequenceFound(ali, data);
						testresults += "  - Sequence name found in alignment: "; 
						testresults += seqID!=-1?"OK\n":"FAILED\n";

						if(seqID!=-1&&dataOK&&dataexists){
							testresults += "  - Data not longer than sequence: "; 
							String newseq = new String(ali.getSequences().get(seqID));
							newseq.replaceAll(gapstring, "");
							testresults += sequenceMatchesProbData(newseq.length(), data)?"OK\n":"FAILED\n";							
						}
					}
				}
				else if(data.getType()==2){
					int seqID = -1;
					boolean dataexists = false; 
					testresults += "- Data identifier: " + data.getiD() + "\n";
					if(data.getFileName()!=null){
						dataexists = attemptToLoadData(data);
						testresults += "  - Data file can be opened: "; 
						testresults += dataexists?"OK\n":"FAILED\n";
					}

					
					if(alignmentOK){
						if(data.getSequenceName()!=null){
							seqID = sequenceFound(ali, data);
							testresults += "  - Sequence name found in alignment: "; 
							testresults += seqID!=-1?"OK\n":"FAILED\n";
						}

						if(seqID!=-1&&dataOK&&dataexists){
							testresults += "  - Data not longer than sequence and no forced pairs are too close: "; 
							String newseq = new String(ali.getSequences().get(seqID));
							newseq.replaceAll(gapstring, "");
							testresults += sequenceMatchesConstraintData(newseq.length(), data)?"OK\n":"FAILED\n";			
						}
						if(data.getContactDistance()>0&&seqID!=-1){
							testresults += "  - Contact distance not longer than alignment: "; 
							testresults += data.getContactDistance()<=ali.getSequences().get(seqID).length()?"OK\n":"FAILED\n";
						}

					}
					
					testresults += "  - NOTE: not checking for pseudoknotted constraints or conflicting data! (Please do that yourself)\n"; 

				}
				else{
					testresults += "- Data identifier: " + data.getiD() + "\n";
					testresults += "  - Unknown data type: " + data.getType() + ", check aborted.";
				}
				testresults += "\n";
			}
		}
		if(PPfoldMain.seqexportname!=null){
			testresults += "\n";
			testresults += "Export check result:\n";
			testresults += "- Sequence name found in alignment: " + (findSequence(ali.getNames(),PPfoldMain.seqexportname)?"OK\n":"FAILED\n");
		}
		
		if(!nothingfailed){
			testresults = SUCCESS_TEXT.concat(testresults);
		}
		else{
			testresults = FAILURE_TEXT.concat(testresults);
		}
				
		return testresults;

	}

	private static long memRequired(Alignment ali) {
		//Calculating that 150 MB is required for 1542 nt
		long size = ali.getSequences().get(0).length();
		long required = (long) ( Math.pow((double)(size)/9158d, 2) * 7270d);
		return required; 
	}

	private static boolean findSequence(List<String> names, String seqexportname) {
			for(int i = 0; i<names.size(); i++){
				if(seqexportname.trim().equals(names.get(i).trim())){
					return true;
				}
			}
			nothingfailed=false;
			return false;
	}

	private static int[] checkProbDataFormat(DataInfo datain) {
		String filename = datain.getFileName();
		BufferedInputStream stream = null;
		int [] readdata_index = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			dataOK=false;
			nothingfailed=false;
			return null;
		}
		
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
				readdata_index = new int[data_size];
				float [] readdata_data1 = new float[data_size];
				float[] readdata_data2 = new float[data_size];
				// Create a pattern to match different kinds of separators
		        Pattern p = Pattern.compile("[,\\s]+");
		        //DO NOT ignore first line... 
				for(int i = 0; i<lines.length; i++){
					line = lines[i];
					String splitline[] = p.split(line.trim());
					if(splitline.length == 3){
						readdata_index[i] = Integer.valueOf(splitline[0])-1; //SHAPE data files are numbered from 1; Java numbers from 0
						readdata_data1[i] = Float.valueOf(splitline[1]);
						readdata_data2[i] = Float.valueOf(splitline[2]);
					}
				}
			}
			catch(Exception e){
				dataOK=false;
				nothingfailed=false;
				return null;
			}
		}
		return readdata_index; 
	}

	private static boolean attemptToLoadData(DataInfo data) {
		String filename = data.getFileName();
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			dataOK=false;
			nothingfailed=false;
			return false;
		}
		return true; 
	}

	private static boolean sequenceMatchesProbData(int length, DataInfo data) {
		// TODO Auto-generated method stub
		return false;
	}

	private static boolean sequenceMatchesConstraintData(int length,
			DataInfo data) {
		ArrayList<int[]> forceArray = new ArrayList<int[]>();
		ArrayList<int[]> prohibitArray = new ArrayList<int[]>();

		//read the data 
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(data.getFileName()));
		} catch (FileNotFoundException e) {
			dataOK=false;
			nothingfailed=false;
			return false;
		}
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
				
				for(int[] cns:forceArray){
					if(cns[0]>length || cns[1]>length || Math.max(cns[1],cns[0]) - Math.min(cns[1], cns[0]) < 4){
						dataOK=false;
						nothingfailed=false;
						return false;
					}
				}
				for(int[] cns:prohibitArray){
					if(cns[0]>length || cns[1]>length){
						dataOK=false;
						nothingfailed=false;
						return false;
					}
				}
				return true; 
			}
			catch(Exception e){
				dataOK=false;
				nothingfailed=false;
				return false;
			}
		}
		
		return false;
	}

	private static boolean sequenceMatchesBarData(int seqlength,
			int[] indices) {
		for(int i = 0; i<indices.length; i++){
			if(indices[i]>seqlength){
				dataOK=false;
				nothingfailed=false;
				return false;
			}
		}
		return true; 
	}

	private static int[] checkBarDataFormat(DataInfo datain) {
		String filename = datain.getFileName();
		int[] indices = null;
		float[] data = null;
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			dataOK=false;
			nothingfailed=false;
			return null;
		}
		
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
				indices = new int[data_size];
				data = new float[data_size];
				// Create a pattern to match different kinds of separators
		        Pattern p = Pattern.compile("[,\\s]+");
				for(int i = 0; i<lines.length; i++){
					line = lines[i];
					String splitline[] = p.split(line.trim());
					if(splitline.length == 2){
						indices[i] = Integer.valueOf(splitline[0])-1; //SHAPE data files are numbered from 1; Java numbers from 0
						data[i] = Float.valueOf(splitline[1]);
					}
					else{
						dataOK=false;
						nothingfailed=false;
						return null;
					}
				}
			}
			catch(Exception e){
				dataOK=false;
				nothingfailed=false;
				return null;
			}
		}
		return indices; 
	}

	private static boolean attemptToLoadBarData(DataInfo data) {
		String filename = data.getFileName();
		BufferedInputStream stream = null;
		try {
			stream = new BufferedInputStream(new FileInputStream(filename));
		} catch (FileNotFoundException e) {
			dataOK=false;
			nothingfailed=false;
			return false;
		}
		return true; 
	}

	private static boolean checkDistProbs(ExtraDataBars dist) {
		float [] distUnpaired = dist.getDistributionUnpaired();
		float [] distPaired = dist.getDistributionPaired();
		float [] limits = dist.getDistributionLimits();
		
		if(!(limits.length==distPaired.length&&limits.length==distUnpaired.length)){
			reason = "Format error";
			System.out.println(limits.length);
			System.out.println(distPaired.length);
			System.out.println(distUnpaired.length);
			distOK=false;
			nothingfailed=false;
			return false;
		}
		
		for(int i = 1; i<limits.length; i++){
			if(limits[i]<=limits[i-1]){
				reason = "Limits are not increasing between " + limits[i-1] + " and " + limits[i];
				distOK=false;
				nothingfailed=false;
				return false;
			}
		}
		
		float sum = 0;
		for(int i = 0; i<limits.length; i++){
			if(distUnpaired[i]>=0){
				sum += distUnpaired[i];
			}
			else{
				reason = "Negative probability " + distPaired[i];
				distOK=false;
				nothingfailed=false;
				return false; 
			}
		}
		if(Math.abs(1f-sum)>0.00001){
			reason = "Unpaired probabilities sum to " + sum;
			distOK=false;
			nothingfailed=false;
			return false;
		}
		
		sum = 0;
		for(int i = 0; i<limits.length; i++){
			if(distPaired[i]>=0){
				sum += distPaired[i];
			}
			else{
				reason = "Negative probability " + distPaired[i];
				distOK=false;
				nothingfailed=false;
				return false; 
			}
		}
		if(Math.abs(1f-sum)>0.00001){
			reason = "Paired probabilities sum to " + sum;
			distOK=false;
			nothingfailed=false;
			return false;
		}
		
		return true;
	}

	private static boolean checkDistHeader(DataInfo data) {
		//read the SHAPE data 
		 BufferedInputStream stream;
			if(data.getDistFileName().equals(PPfoldMain.defaultDataDistfile)){
				stream = new BufferedInputStream(Thread.currentThread().
						getContextClassLoader().getResourceAsStream(PPfoldMain.defaultDataDistfile));
			}
			else{
				try {
					stream = new BufferedInputStream(new FileInputStream(data.getDistFileName()));
				} catch (FileNotFoundException e) {
					nothingfailed = false; 
					distOK = false;
					return false;				}
			}
			try{
					int l = stream.available();
					byte[] bytes = new byte[l];
					stream.read(bytes);
					stream.close();
					String data_string = new String(bytes);
					String[] lines = data_string.split("\n");
					String line = lines[0]; //First line
					String[] splitline = line.split("\\s+");
					if(!( splitline[0].startsWith("lower_bound") && 
							splitline[1].startsWith("P_density_paired") && 
							splitline[2].startsWith("P_density_unpaired"))){
						nothingfailed = false; 
						distOK = false;
						return false;
					}
					else{return true;}
			}
			catch(Exception e){
				nothingfailed = false; 
				distOK = false;
				return false;
			}
	}

	private static int sequenceFound(Alignment ali, DataInfo data) {
		Integer seqID = -1;
		for(int i = 0; i<ali.getNames().size(); i++){
			if(data.getSequenceName().trim().equals(ali.getNames().get(i).trim())){
				seqID = i;
				break;
			}
		}
		if(seqID==-1){dataOK=false; nothingfailed=false;return seqID;}
		else{return seqID;}
	}

	private static boolean attemptToAccessOutput() {
		if(PPfoldMain.outputdir==null){
			nothingfailed=false;
			return false;
		}
		File file=new File(PPfoldMain.outputdir);
		if(!file.exists()){nothingfailed=false;}
		return file.exists();
	}

	private static boolean attemptToWriteOutput() {
		File sample = null;
		if(PPfoldMain.outputdir!=null){
			sample = new File(PPfoldMain.outputdir,"tmp"); 
		}
		else{nothingfailed=false;return false;}
		try{
		      sample.createNewFile();
		      sample.delete();
		      return true;
		}
		catch(IOException e)
		{
			nothingfailed=false;
		      return false;
		}
	}

	private static String namesMatchTreeAlignment(Alignment ali, Tree tree) {
		for (int i = 0; i < ali.getNames().size(); i++) {
			// finds the node corresponding to the rownumber
			// rownumber = corresponds to sequences.
			Node node = tree.findSlowlyNodeWithName(ali.getNames().get(i));
			if (node == null) {
				nothingfailed=false;
				return ali.getNames().get(i);
			}
		}
		nothingfailed=false;
		return null;
	}

	private static Tree attemptToLoadTree() {
		 try{
			 File file=new File(PPfoldMain.treefilename);
				if(!file.exists()){nothingfailed=false;return null;}
				else if(file.length()>1048576){
					nothingfailed=false;
					fileTooLong=true;
					return null;
				}
	        	Tree tree = NewickReader.readNewick(PPfoldMain.treefilename);
	        	return tree;
		 }
		 catch(Exception e){nothingfailed=false;return null;}
	}

	private static int lengthOK(Alignment ali) {
		try{
		int size = ali.getSequences().get(0).length();
		for(int i = 0; i<ali.getSequences().size();i++){
			if(ali.getSequences().get(i).length()!=size){
				alignmentOK=false;nothingfailed=false;
				return i;
			}
		}
		return 0;
		}
		catch(Exception e){
			return -1;
		}
	}

	private static boolean namesMatchSequences(Alignment ali) {
		boolean val = (ali.getNames().size()==ali.getSequences().size())&&(ali.getNames().size()!=0);
		if(!val){alignmentOK=false;nothingfailed=false;}
		return val;
	}

	private static String noIllegalSymbols(Alignment ali) {
		for(int i = 0; i<ali.getSequences().size();i++){
			for(int j=0; j<ali.getSequences().get(i).length(); j++){
				char thischar = ali.getSequences().get(i).charAt(j);
				thischar = Character.toLowerCase(thischar);
				if (thischar != 'a' && thischar != 'u' && thischar != 't'
						&& thischar != 'g' && thischar != 'c'
						&& thischar != 'r' && thischar != 'y'
						&& thischar != 's' && thischar != 'w'
						&& thischar != 'k' && thischar != 'm'
						&& thischar != 'b' && thischar != 'd'
						&& thischar != 'h' && thischar != 'v'
								&& thischar != 'n'
						&& !MatrixTools.isGap(thischar)){
					alignmentOK=false;nothingfailed=false;
					return String.valueOf(thischar);
				}
			}
		}
		return null;
	}

	private static Alignment attemptToLoadAlignment() {
		 try{
			 	File file=new File(PPfoldMain.alignmentfilename);
				if(!file.exists()){nothingfailed=false;alignmentOK=false;return null;}
				else if(file.length()>1048576){
					alignmentOK=false;
					nothingfailed=false;
					fileTooLong=true;
					return null;
				}
			 
	        	Alignment align = AlignmentReader.readAlignment(PPfoldMain.alignmentfilename);
	        	return align;
		 }
		 catch(Exception e){alignmentOK=false;nothingfailed=false;return null;}
	}
	
	private static ExtraDataBars attemptToLoadDataDist(DataInfo data) {
		 try{					
			 BufferedInputStream shapeDistReader;
			if(data.getDistFileName().equals(PPfoldMain.defaultDataDistfile)){
				shapeDistReader = new BufferedInputStream(Thread.currentThread().
						getContextClassLoader().getResourceAsStream(PPfoldMain.defaultDataDistfile));
			}
			else{
				shapeDistReader = new BufferedInputStream(new FileInputStream(data.getDistFileName()));
			}
			ExtraDataBars sequenceData = ExtraDataBars.readDistTable_toStream(shapeDistReader);
	        return sequenceData;
		 }
		 catch(Exception e){distOK=false;nothingfailed=false;return null;}
	}
	
	

}
