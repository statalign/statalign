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


public class RNAalifold {
	
	public static String executable = "/home/michael/Downloads/ViennaRNA-2.0.7/Progs/RNAalifold";
	
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
				res = RNAalifold.fold(sequences, sequenceNames,newparams);
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
				res = RNAalifold.fold(sequences, sequenceNames,oldparams);
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
	
	public static RNAalifoldResult fold(List<String> sequences, List<String> sequenceNames, String arguments) throws Exception
	{
		return fold(sequences, sequenceNames, arguments, true);
	}
	
	public static RNAalifoldResult fold(List<String> sequences, List<String> sequenceNames, String arguments, boolean useMatrix) throws Exception
	{
		if(useOldParams)
		{
			arguments.replaceAll(" --cfactor ", " -cv ");
			arguments.replaceAll(" --nfactor ", " -nc ");
			//System.out.println(arguments);
		}
		
		try
		{
			System.out.println(arguments);
			File tempClustalFile = new File("temp.clustalw");
			saveClustalW(sequences, sequenceNames, tempClustalFile);
			String args = executable + " " + "-p --bppmThreshold=0 "+arguments+" "+tempClustalFile.getAbsolutePath();
			if(useOldParams)
			{
				args = executable + " " + "-p "+arguments+" "+tempClustalFile.getAbsolutePath();
			}
			Process p = Runtime.getRuntime().exec(args);
			InputStream is = p.getErrorStream();
			BufferedReader buffer = new BufferedReader(new InputStreamReader(is));
			String textline = null;
			String errorString = "";
			boolean first = true;
			while ((textline = buffer.readLine()) != null) {
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
			int exitCode = p.waitFor();
			//System.out.println(exitCode);
			if(exitCode != 0)
			{
				System.err.println("RNAalifold generated the following error during execution:" +
						"\n\"" + errorString+"\"");
				return null;
				//throw new Exception();
				//System.out.println("The following error occured:");
				//System.err.print(errorString);
			}

			RNAalifoldResult result = new RNAalifoldResult();
			if(useMatrix)
			{
				result.matrix = loadBasePairProbMatrix(new File("alidot.ps"), sequences.get(0).length());
			}
			result.pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketString(loadDotBracketStructure(new File("alifold.out")));
			return result;
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
		}
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