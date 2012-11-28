package statalign.postprocess.utils;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import statalign.postprocess.plugins.RNAalifoldResult;


/**
 * Given an alignment of sequence this class shells to the RNAalifold
 * appropriately and returns a RNAalifold secondary structure prediction.
 * 
 */
public class RNAalifold {
	
	public static String executable = "RNAalifold.exe";
	
	static boolean useOldParams = false;
	
	public static boolean checkRNAalifold()
	{
		try
		{
			ArrayList<String> sequences = new ArrayList<String>();
			ArrayList<String> sequenceNames = new ArrayList<String>();
			sequences.add("GGGUGCUUGAAGCUGUCUGCUUUAAGUGCUUGCA----UCAGGCUGAGAGUAGGCAGAGAAAAGCCCCGUAUCA-----A----------------UGUUAAUCAAUACGAGGC-CCUCUGUAAUG");
			sequences.add("GGGUGCUUGAGGCUGUCUGCCUCGGG------CAUGCC---ACCGUAAGGCAGACAGAGAAAAGCCCCAGUUAACAUUACGCGUCCUGCAAGACGCCUAACAUUAAUCUGAGGC-CAAUUU-CAUG");
			sequenceNames.add("a");
			sequenceNames.add("b");
			
			useOldParams = false;
			String newparams = " -T " + 37 +" --cfactor " +  1 + " --nfactor " + 1 + " ";
			RNAalifoldResult res = null;
			try
			{
				res = RNAalifold.fold(sequences, sequenceNames,newparams, true);
			}
			catch(Exception ex)
			{
				System.err.println("The following error occured with RNAalifold: " + ex.getMessage());
			}
			//System.out.println("HERE " + res);
			if(res != null)
			{
				return true;
			}
			else
			{
				String oldparams = " -T " + 37 +" -cv " +  1 + " -nc " + 1 + " ";
				useOldParams = true;
				res = RNAalifold.fold(sequences, sequenceNames,oldparams, true);
				return res != null;
			}
		}
		catch(Exception ex)
		{
		System.err.println("The following error occured with RNAalifold: " + ex.getMessage());
		}

		useOldParams = false;
		
		return false;
	}
	

	public static RNAalifoldResult fold(List<String> sequences, List<String> sequenceNames, String arguments)  throws Exception
	{
		return fold(sequences, sequenceNames, arguments, true, false);
	}
	
	public static RNAalifoldResult fold(List<String> sequences, List<String> sequenceNames, String arguments, boolean noErrorMessages) throws Exception
	{
		return fold(sequences, sequenceNames, arguments, true, noErrorMessages);
	}
	
	public static RNAalifoldResult fold(List<String> sequences, List<String> sequenceNames, String arguments, boolean useMatrix, boolean noErrorMessages) throws Exception
	{
		if(useOldParams)
		{
			//System.out.println("using old");
			//System.out.println("C"+arguments);
			arguments = arguments.replaceAll(" --cfactor ", " -cv ");
			arguments = arguments.replaceAll(" --nfactor ", " -nc ");
			//System.out.println("D"+arguments);
		}
		
		String tempPath = System.getProperty("java.io.tmpdir")+"/";		
		File tempClustalFile = new File(tempPath+"temp.clustalw");
		saveClustalW(sequences,sequenceNames,tempClustalFile);
		try
		{	
			
			System.out.println(arguments);
			
			System.out.println(tempClustalFile);
			String args = executable + " " + "-p "+arguments;
			String file = tempClustalFile.getAbsolutePath();
			if(useOldParams)
			{
				args = executable + " " + "-p "+arguments;
			}
			
			
			// create a process builder to execute in temporary directory
			ProcessBuilder processBuilder = new ProcessBuilder();			
			processBuilder.directory(new File(tempPath));
			ArrayList<String> commands = new ArrayList<String>();
			String [] split = args.split("(\\s)+");
			for(int i = 0 ; i < split.length ; i++)
			{
				commands.add(split[i]);
			}
			commands.add(file);
			//System.out.println(commands);
			processBuilder.command(commands);
			
			Process p = processBuilder.start();
			
			InputStream is = p.getErrorStream();
			BufferedReader buffer = new BufferedReader(new InputStreamReader(is));
			String textline = null;
			String errorString = "";
			boolean first = true;
			while ((textline = buffer.readLine()) != null) {

				System.out.println(textline);
				if(first)
				{
					first = false;
				}
				else
				{
					errorString += "\n";
				}
				errorString += textline;
			}
			buffer.close();
			//System.out.println("A");
			int exitCode = p.waitFor();
			//System.out.println(exitCode);
			if(exitCode != 0)
			{
				if(!noErrorMessages)
				{
				System.err.println("RNAalifold generated the following error during execution:" +
						"\n\"" + errorString+"\"");
				}
				return null;
				//throw new Exception();
				//System.out.println("The following error occured:");
				//System.err.print(errorString);
			}
			//System.out.println("B");
			RNAalifoldResult result = new RNAalifoldResult();
			if(useMatrix)
			{
				File matrixFile = new File(tempPath+"alidot.ps");
				result.matrix = loadBasePairProbMatrix(matrixFile, sequences.get(0).length());
				matrixFile.delete();
			}
			File outFile = new File(tempPath+"alifold.out");
			result.pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketString(loadDotBracketStructure(outFile));
			outFile.delete();
			tempClustalFile.delete();
			return result;
		}
		catch(Exception ex)
		{			
			//ex.printStackTrace();
		}
		tempClustalFile.delete();
		return null;
	}
	
	public static void saveClustalW(List<String> sequences, List<String> sequenceNames, File outFile)
	{
		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile));
			buffer.write("CLUSTAL W(1.60) multiple sequence alignment\n");
			buffer.write("\n");
			for(int i = 0 ; i < sequences.size() ; i++)
			{
				buffer.write(sequenceNames.get(i)+"\t"+sequences.get(i)+" "+sequences.get(i).length()+"\n");
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	public static String loadDotBracketStructure(File alifoldOut)
	{
		String dbs = null;
		
		try
		{
			BufferedReader buffer = new BufferedReader(new FileReader(alifoldOut));
			String textline = null;
			
			while((textline = buffer.readLine()) != null)
			{
				dbs = textline;
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
		return dbs;
	}
	
	public static double[][] loadBasePairProbMatrix(File basePairFile, int length)
	{
		double [][] matrix = new double[length][length];
		try
		{
			BufferedReader buffer = new BufferedReader(new FileReader(basePairFile));
			String textline = null;
			boolean cont = false;
			while((textline = buffer.readLine()) != null)
			{
				if(cont)
				{
					String [] split = textline.split("(\\s)+");
					if(split.length >= 7 && split[6].endsWith("ubox"))
					{
						int x = Integer.parseInt(split[3]) - 1;
						int y = Integer.parseInt(split[4]) - 1;
						double prob = Double.parseDouble(split[5]);
						matrix[x][y] = prob;
						matrix[y][x] = prob;
						//System.out.println(x+"\t"+y+"\t"+prob);
					}
				}
				else
				if(textline.startsWith("drawgrid"))
				{
					cont = true;
				}
			}		
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
		
		return matrix;
	}
	
	public static void main(String[] args)
	{
		try
		{
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> sequenceNames = new ArrayList<String>();
		RNAFoldingTools.loadFastaSequences(new File("/home/michael/Dropbox/RNA and StatAlign/Distance/Datasets2/TestRNAData1_5seqs.dat.fas"), sequences, sequenceNames);
		RNAalifold.executable="/home/michael/Downloads/ViennaRNA-2.0.7/Progs/RNAalifold";
		RNAalifoldResult res = RNAalifold.fold(sequences, sequenceNames,"-T 10");
		System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(res.pairedSites));
		
		res = RNAalifold.fold(sequences, sequenceNames,"-T 60");
		System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(res.pairedSites));
		//RNAalifold.saveClustalW(sequences, sequenceNames, new File("/home/michael/Desktop/temp.clustalw"));
		}
		catch(Exception ex)
		{
			
		}
	}
}