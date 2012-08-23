package statalign.postprocess.plugins.benchmarks;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import com.ppfold.algo.FuzzyAlignment;
import com.ppfold.algo.ResultBundle;

import statalign.postprocess.plugins.PPFold;
import statalign.postprocess.plugins.contree.hash.HashTable;
import statalign.postprocess.utils.Mapping;
import statalign.postprocess.utils.RNAFoldingTools;


public class Benchmarks 
{
	public static void main(String[] args) {
		Benchmarks.automatedTests();
		//Benchmarks.testVariation2();
		//new Benchmarks().performEntropy();
		//Benchmarks.automatedTest2();
		System.exit(0);
		/*
		//Benchmarks.testData();
		Benchmarks.automatedTest2();
		//.automatedTest();
		System.exit(0);
		*/
		
		String dir = "/home/michael/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		File [] files = new File(dir).listFiles();
		for(int i = 0 ; i < files.length ; i++)
		{
			if(files[i].getName().toLowerCase().endsWith(".dat"))
			{
				ExperimentalData expData = Benchmarks.loadExperimentalStructure(files[i]);
				saveAsFasta(expData, new File(dir+files[i].getName()+".fas"));
			}
		}
		
		dir = "/home/michael/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		files = new File(dir).listFiles();
		for(int i = 0 ; i < files.length ; i++)
		{
			if(files[i].getName().toLowerCase().endsWith(".ct"))
			{
				try
				{
					BufferedWriter buffer = new BufferedWriter(new FileWriter(new File(dir+files[i].getName()+".dbn")));
					buffer.write(RNAFoldingTools.getDotBracketStringFromCtFile(files[i])+"\n");
					buffer.close();
				}
				catch(IOException ex)
				{
					ex.printStackTrace();
				}
			}
		}
		
		/*
		ExperimentalData expData = Benchmarks.loadExperimentalStructure(new File("C:/Oxford/TestRNAData.tar/TestRNAData/TestRNAData1.dat"));
		ExperimentalData predictedData = Benchmarks.loadExperimentalStructure(new File("C:/Oxford/TestRNAData.tar/TestRNAData/TestRNAData1Predicted.dat"));
		System.out.println(calculateSensitivity(expData.pairedSites, predictedData.pairedSites));
		System.out.println(calculatePPV(expData.pairedSites, predictedData.pairedSites));
		System.out.println(calculateFScore(expData.pairedSites, predictedData.pairedSites));
		*/
		
	}
	
	public static void testData()
	{
		String name = "TestRNAData1";

		String dir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/";
		String resultsDir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/Results0/";
		File experimentalFile = new File(dir+name+".dat");
		
		
		ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
		StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res"));
		
		//System.out.println(statalignResult.sequence+"");
		int [] pairedSitesExperimental = projectPairedSites(statalignResult.sequence, experimentalData.pairedSites);
		int [] pairedSitesStatAlign = statalignResult.pairedSites;
		File ppfoldData = new File("/home/michael/Dropbox/RNA and StatAlign/TestRNAData/StatAlign/1.ct");
		
		/*
		double [][] bpMatrix = RNAFoldingTools.loadMatrix(new File("/home/michael/Dropbox/RNA and StatAlign/TestRNAData/StatAlign/1.bp"));
		float [][] floatMatrix = new float[bpMatrix.length][bpMatrix[0].length];
		for(int i = 0; i < bpMatrix.length ; i++)
		{
			for(int j = 0; j < bpMatrix[0].length ; j++)
			{
				floatMatrix[i][j] = (float) bpMatrix[i][j];
			}
		}
		String s = "GGGCGCCCGAGGCCGCCCGCCCCGGGCACGCCACCGCAAG------GCAGACAGAGAAAAGCCCCAGCCAACACCACGCGCCCCGCAAGACGCCCAACACCAA-CCCGAGGCCCAAC-CCACGCCCCACAAACGCAGGCCAGCCCCCCACGCGCCGAAAGGCAAG---GAGAAGCAGGCCACGAAG";
		float [][] projectMatrix = Mapping.projectMatrix(s, floatMatrix, '-');
		double [][] doubleMatrix = new double[projectMatrix.length][projectMatrix[0].length];
		for(int i = 0; i < doubleMatrix.length ; i++)
		{
			for(int j = 0; j < doubleMatrix[0].length ; j++)
			{
				doubleMatrix[i][j] = floatMatrix[i][j];
			}
		}
		*/
		//int [] pairedSitesPPfold = RNAFoldingTools.getPosteriorDecodingConsensusStructure(doubleMatrix);
		
		
		System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData)));
		int [] pairedSitesPPfold = projectPairedSites("GGGCGCCCGAGGCCGCCCGCCCCGGGCACGCCACCGCAAG------GCAGACAGAGAAAAGCCCCAGCCAACACCACGCGCCCCGCAAGACGCCCAACACCAA-CCCGAGGCCCAAC-CCACGCCCCACAAACGCAGGCCAGCCCCCCACGCGCCGAAAGGCAAG---GAGAAGCAGGCCACGAAG", RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));

		//File ppfoldData = new File("/host/Oxford/test1_mpd.ct");
		//int [] pairedSitesPPfold = projectPairedSites("GGGCGCCCGAGGCCGCCCGCCCCGGGCACGCCACCGCAAG------GCAGACAGAGAAAAGCCCCAGCCAACACCACGCGCCCCGCAAGACGCCCAACACCAA-CCCGAGGCCCAAC-CCACGCCCCACAAACGCAGGCCAGCCCCCCACGCGCCGAAAGGCAAG---GAGAAGCAGGCCACGAAG", RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
		//printPairs(pairedSitesPPfold);
		System.out.println("LENGTH:"+pairedSitesStatAlign.length);
		
		System.out.println(">"+name + " (" + (statalignResult.sequence.replaceAll("-", "").length())+")");
		System.out.println("  "+statalignResult.sequence.replaceAll("-", ""));
		System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
		System.out.println("S:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesStatAlign));
		System.out.println("  "+statalignResult.sequence.replaceAll("-", ""));
		System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
		System.out.println("P:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesPPfold));
		
		
		double sensExpStat = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesStatAlign);
		double sensExpPPfold = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesPPfold);
		double ppvExpStat = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesStatAlign);
		double ppvExpPPfold = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesPPfold);
		double fscExpStat = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlign);
		double fscExpPPfold = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesPPfold);
		
		System.out.println(RNAFoldingTools.pad("", 14)+RNAFoldingTools.pad("StatAl", 10)+RNAFoldingTools.pad("PPfold", 10));
		System.out.println(RNAFoldingTools.pad("Senstivity", 14)+RNAFoldingTools.pad(sensExpStat+"", 6)+"    "+RNAFoldingTools.pad(sensExpPPfold+"", 6));
		System.out.println(RNAFoldingTools.pad("PPV", 14)+RNAFoldingTools.pad(ppvExpStat+"", 6)+"    "+RNAFoldingTools.pad(ppvExpPPfold+"", 6));
		System.out.println(RNAFoldingTools.pad("F-score", 14)+RNAFoldingTools.pad(fscExpStat+"", 6)+"    "+RNAFoldingTools.pad(fscExpPPfold+"", 6));
		System.out.print("StatAl: ");
		printValues(pairedSitesExperimental, pairedSitesStatAlign);
		System.out.print("PPfold: ");
		printValues(pairedSitesExperimental, pairedSitesPPfold);
		System.out.println("---------------------------------------------------------------------");
		System.out.println();
	}
	
	public static void automatedTest()
	{
		String dir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/";
		String resultsDir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/Results0/";
		File [] files = new File(resultsDir).listFiles();
		for(int i = 0 ; i < files.length ; i++)
		{
			String fullName = files[i].getName();
			if(fullName.endsWith(".dat.res"))
			{
				String name = fullName.substring(0, fullName.length()-8);

				
				//File files = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/").listFiles();
				File experimentalFile = new File(dir+name+".dat");
				//File ourData = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNADATA1OURS");
				//File statalignResultFile = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNAData1.dat.txt");
				File ppfoldData = new File(dir+name+".dat.ct");
				
				ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
				StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res"));
				
				
				String mappingSeq = "";
				for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
				{
					if(experimentalData.sequences.get(j).replaceAll("-", "").equals(statalignResult.sequence.replaceAll("-", "")))
					{
						mappingSeq = experimentalData.sequences.get(j);
					}
				}
				
				//System.out.println(mappingSeq);
				//System.out.println(statalignResult.sequence);
				
				int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
				int [] pairedSitesStatAlign = statalignResult.pairedSites;
				
				int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
				//printPairs(pairedSitesExperimental);
				//printPairs(pairedSitesStatAlign);
				//printPairs(pairedSitesPPfold);
				

				//System.out.println("X:"+statalignResult.sequence);
				//System.out.println("X:"+experimentalData.sequences.get(2));
				System.out.println(">"+name + " (" + (statalignResult.sequence.replaceAll("-", "").length())+")");
				System.out.println("  "+statalignResult.sequence.replaceAll("-", ""));
				System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
				System.out.println("S:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesStatAlign));
				System.out.println("  "+statalignResult.sequence.replaceAll("-", ""));
				System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
				System.out.println("P:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesPPfold));
				
				
				double sensExpStat = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesStatAlign);
				double sensExpPPfold = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesPPfold);
				double ppvExpStat = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesStatAlign);
				double ppvExpPPfold = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesPPfold);
				double fscExpStat = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlign);
				double fscExpPPfold = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesPPfold);
				
				System.out.println(RNAFoldingTools.pad("", 14)+RNAFoldingTools.pad("StatAl", 10)+RNAFoldingTools.pad("PPfold", 10));
				System.out.println(RNAFoldingTools.pad("Senstivity", 14)+RNAFoldingTools.pad(sensExpStat+"", 6)+"    "+RNAFoldingTools.pad(sensExpPPfold+"", 6));
				System.out.println(RNAFoldingTools.pad("PPV", 14)+RNAFoldingTools.pad(ppvExpStat+"", 6)+"    "+RNAFoldingTools.pad(ppvExpPPfold+"", 6));
				System.out.println(RNAFoldingTools.pad("F-score", 14)+RNAFoldingTools.pad(fscExpStat+"", 6)+"    "+RNAFoldingTools.pad(fscExpPPfold+"", 6));
				System.out.print("StatAl: ");
				printValues(pairedSitesExperimental, pairedSitesStatAlign);
				System.out.print("PPfold: ");
				printValues(pairedSitesExperimental, pairedSitesPPfold);
				System.out.println("---------------------------------------------------------------------");
				System.out.println();
			}
		}
	}
	
	/*public static void automatedTest2()
	{
		String dir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/";
		String resultsDir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/Results3/";
		File [] files = new File(resultsDir).listFiles();
		for(int i = 0 ; i < files.length ; i++)
		{
			String fullName = files[i].getName();
			if(fullName.endsWith(".dat.res"))
			{
				String name = fullName.substring(0, fullName.length()-8);

				
				//File files = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/").listFiles();
				File experimentalFile = new File(dir+name+".dat");
				//File ourData = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNADATA1OURS");
				//File statalignResultFile = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNAData1.dat.txt");
				File ppfoldData = new File(dir+name+".dat.ct");
				
				ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
				StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res"));
				StatAlignResult statalignWeightedResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res.weighted"));
				System.out.println(name);
				StatAlignResult mpdResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res.mpd"));
				
				String mappingSeq = "";
				for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
				{
					if(experimentalData.sequences.get(j).replaceAll("-", "").equals(statalignResult.sequence.replaceAll("-", "")))
					{
						mappingSeq = experimentalData.sequences.get(j);
					}
				}
				
				//System.out.println(mappingSeq);
				//System.out.println(statalignResult.sequence);
				
				int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
				int [] pairedSitesStatAlign = statalignResult.pairedSites;
				int [] pairedSitesStatAlignWeighted = statalignWeightedResult.pairedSites;
				int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
				int [] pairedSitesMPD = mpdResult.pairedSites;
				//printPairs(pairedSitesExperimental);
				//printPairs(pairedSitesStatAlign);
				//printPairs(pairedSitesPPfold);
				

				//System.out.println("X:"+statalignResult.sequence);
				//System.out.println("X:"+experimentalData.sequences.get(2));
				System.out.println(">"+name + " (" + (statalignResult.sequence.replaceAll("-", "").length())+")");
				System.out.println("  "+statalignResult.sequence.replaceAll("-", ""));
				System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
				System.out.println("S:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesStatAlign));
				System.out.println("  "+statalignResult.sequence.replaceAll("-", ""));
				System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
				System.out.println("P:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesPPfold));
				
				
				double sensExpStat = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesStatAlign);
				double sensExpPPfold = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesPPfold);
				double ppvExpStat = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesStatAlign);
				double ppvExpPPfold = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesPPfold);
				double fscExpStat = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlign);
				double fscExpPPfold = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesPPfold);
				double fscExpStatWeighted =Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlignWeighted);
				double fscExpStatMPD=Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesMPD);
				
				System.out.println(RNAFoldingTools.pad("", 14)+RNAFoldingTools.pad("StatAl", 10)+RNAFoldingTools.pad("PPfold", 10));
				System.out.println(RNAFoldingTools.pad("Senstivity", 14)+RNAFoldingTools.pad(sensExpStat+"", 6)+"    "+RNAFoldingTools.pad(sensExpPPfold+"", 6));
				System.out.println(RNAFoldingTools.pad("PPV", 14)+RNAFoldingTools.pad(ppvExpStat+"", 6)+"    "+RNAFoldingTools.pad(ppvExpPPfold+"", 6));
				System.out.println(RNAFoldingTools.pad("F-score", 14)+RNAFoldingTools.pad(fscExpStat+"", 6)+"    "+RNAFoldingTools.pad(fscExpPPfold+"", 6));
				System.out.print("StatAl: ");
				printValues(pairedSitesExperimental, pairedSitesStatAlign);
				System.out.print("PPfold: ");
				printValues(pairedSitesExperimental, pairedSitesPPfold);
				System.out.println("---------------------------------------------------------------------");
				System.out.println();
				
				//System.out.println("XXXXXXXXXXXXXX"+resultsDir+name+".folds");
				if(new File(resultsDir+name+".folds").exists())
				{
					ArrayList<String> structures = PPFold.loadFolds(new File(resultsDir+name+".folds"), 4);
					ArrayList<String> values = new ArrayList<String>();
					for(int k= 0 ; k < structures.size() ; k++)
					{
						String val = "" + Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(structures.get(k)));
						values.add(val);
					}
					//System.out.println(values);
					try
					{
						BufferedWriter buffer = new BufferedWriter(new FileWriter(resultsDir+name+".hist"));
						buffer.write("ST="+fscExpStat+"\n");
						buffer.write("STW="+fscExpStatWeighted+"\n");
						buffer.write("MPD="+fscExpStatMPD+"\n");
						buffer.write("PP="+fscExpPPfold+"\n");						
						for(int l = 0 ; l < values.size() ; l++)
						{
							double val = Double.parseDouble(values.get(l));
							if(Double.isNaN(val))
							{
								val = 0;
							}
							buffer.write(val+"\n");
						}
						buffer.close();
					}
					catch(IOException ex)
					{
						ex.printStackTrace();
					}
					
					String dir2 = "/home/michael/workspace/StatAlign/";
					System.out.println(new File(dir2+name+".dat.fas.folds_e_obs"));
					if(new File(dir2+name+".dat.fas.folds_e_obs").exists())
					{
						String dbn = PPFold.loadFolds(new File(dir2+name+".dat.fas.folds_e_obs"), 4).get(0);
						double fscSamplingObs = Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(dbn));
						System.out.println("Alignment sampling obs FSC"+fscSamplingObs);
						
						dbn = PPFold.loadFolds(new File(dir2+name+".dat.fas.folds_e_exp"), 4).get(0);
						double fscSamplingExp = Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(dbn));
						System.out.println("Alignment sampling exp FSC"+fscSamplingExp);
						
					}
					//System.out.println(structures);
				}
			}
		}
	}*/
	
	public static void automatedTest2()
	{
		//String dir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		String dir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/";
		//String resultsDir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/Results3/";
		String resultsDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq2/";
		File [] files = new File(resultsDir).listFiles();
		for(int i = 0 ; i < files.length ; i++)
		{
			String fullName = files[i].getName();
			if(fullName.endsWith(".dat.res"))
			{
				String name = fullName.substring(0, fullName.length()-8);
				String smallname = fullName.substring(0, fullName.length()-16);
				System.out.println(name);
				
				String truncname = smallname.replaceAll("_5seqs", "");
				//File files = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/").listFiles();
				File experimentalFile = new File(dir+truncname+".dat");
				//File ourData = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNADATA1OURS");
				//File statalignResultFile = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNAData1.dat.txt");
				File ppfoldData = new File(dir+truncname+".dat.ct");
				
				ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
				StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res"));
				StatAlignResult statalignWeightedResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res.weighted"));
				System.out.println(name);
				StatAlignResult mpdResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res.mpd"));
				
				String mappingSeq = "";
				for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
				{
					if(experimentalData.sequences.get(j).replaceAll("-", "").equals(statalignResult.sequence.replaceAll("-", "")))
					{
						mappingSeq = experimentalData.sequences.get(j);
					}
				}
				
				//System.out.println(mappingSeq);
				//System.out.println(statalignResult.sequence);
				
				int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
				int [] pairedSitesStatAlign = statalignResult.pairedSites;
				int [] pairedSitesStatAlignWeighted = statalignWeightedResult.pairedSites;
				int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
				int [] pairedSitesMPD = mpdResult.pairedSites;
				//printPairs(pairedSitesExperimental);
				//printPairs(pairedSitesStatAlign);
				//printPairs(pairedSitesPPfold);
				

				//System.out.println("X:"+statalignResult.sequence);
				//System.out.println("X:"+experimentalData.sequences.get(2));
				System.out.println(">"+name + " (" + (statalignResult.sequence.replaceAll("-", "").length())+")");
				System.out.println("  "+statalignResult.sequence.replaceAll("-", ""));
				System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
				System.out.println("S:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesStatAlign));
				System.out.println("  "+statalignResult.sequence.replaceAll("-", ""));
				System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
				System.out.println("P:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesPPfold));
				
				
				double sensExpStat = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesStatAlign);
				double sensExpPPfold = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesPPfold);
				double ppvExpStat = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesStatAlign);
				double ppvExpPPfold = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesPPfold);
				double fscExpStat = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlign);
				double fscExpPPfold = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesPPfold);
				double fscExpStatWeighted =Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlignWeighted);
				double fscExpStatMPD=Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesMPD);
				double fscSamplingObs = -1;
				
				System.out.println(RNAFoldingTools.pad("", 14)+RNAFoldingTools.pad("StatAl", 10)+RNAFoldingTools.pad("PPfold", 10));
				System.out.println(RNAFoldingTools.pad("Senstivity", 14)+RNAFoldingTools.pad(sensExpStat+"", 6)+"    "+RNAFoldingTools.pad(sensExpPPfold+"", 6));
				System.out.println(RNAFoldingTools.pad("PPV", 14)+RNAFoldingTools.pad(ppvExpStat+"", 6)+"    "+RNAFoldingTools.pad(ppvExpPPfold+"", 6));
				System.out.println(RNAFoldingTools.pad("F-score", 14)+RNAFoldingTools.pad(fscExpStat+"", 6)+"    "+RNAFoldingTools.pad(fscExpPPfold+"", 6));
				System.out.print("StatAl: ");
				printValues(pairedSitesExperimental, pairedSitesStatAlign);
				System.out.print("PPfold: ");
				printValues(pairedSitesExperimental, pairedSitesPPfold);
				System.out.println("---------------------------------------------------------------------");
				System.out.println();
				
				//String dir2 = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq/";
				String dir2 = resultsDir;
				//System.out.println(new File(dir2+smallname+".dat.fas.folds_e_obs").exists());
				if(new File(dir2+smallname+".dat.fas.folds_e_obs").exists())
				{
					String dbn = PPFold.loadFolds(new File(dir2+smallname+".dat.fas.folds_e_obs"), 4).get(0);
					fscSamplingObs = Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(dbn));
					System.out.println("Alignment sampling obs FSC"+fscSamplingObs);
					
					dbn = PPFold.loadFolds(new File(dir2+smallname+".dat.fas.folds_e_exp"), 4).get(0);
					double fscSamplingExp = Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(dbn));
					System.out.println("Alignment sampling exp FSC"+fscSamplingExp);
					
				}
				
				System.out.println(new File(resultsDir+name+".folds").exists()+"\t"+fscSamplingObs);
				//System.out.println("XXXXXXXXXXXXXX"+resultsDir+name+".folds");
				if(new File(resultsDir+name+".folds").exists() && fscSamplingObs != -1)
				{
					ArrayList<String> structures = PPFold.loadFolds(new File(resultsDir+name+". >folds"), 4);
					ArrayList<String> values = new ArrayList<String>();
					for(int k= 0 ; k < structures.size() ; k++)
					{
						String val = "" + Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(structures.get(k)));
						values.add(val);
					}
					//System.out.println(values);
					try
					{
						BufferedWriter buffer = new BufferedWriter(new FileWriter(resultsDir+name+".hist2"));
						buffer.write("ST="+fscExpStat+"\n");
						buffer.write("STW="+fscExpStatWeighted+"\n");
						buffer.write("MPD="+fscExpStatMPD+"\n");
						buffer.write("PP="+fscExpPPfold+"\n");
						buffer.write("STE="+fscSamplingObs+"\n");
						for(int l = 0 ; l < values.size() ; l++)
						{
							double val = Double.parseDouble(values.get(l));
							if(Double.isNaN(val))
							{
								val = 0;
							}
							buffer.write(val+"\n");
						}
						buffer.close();
					}
					catch(IOException ex)
					{
						ex.printStackTrace();
					}
					
					
					//System.out.println(structures);
				}
			}
		}
	}
	
	public static void performDistanceBenchmarks(Dataset dataset)
	{
		//File distanceFile = new File(System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/9seq2/dist_scores.txt");
		/*
		String dataDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		String resultsDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/9seq2/";
		
					
					File experimentalFile = new File(dataDir+name+".dat");
					File ppfoldData = new File(dataDir+name+".dat.ct");
					
					ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
					
					StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res")); // 
					StatAlignResult statalignWeightedResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res.weighted"));
					StatAlignResult mpdResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res.mpd"));
					
					String mappingSeq = "";
					for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
					{
						if(experimentalData.sequences.get(j).replaceAll("-", "").equals(statalignResult.sequence.replaceAll("-", "")))
						{
							mappingSeq = experimentalData.sequences.get(j);
						}
					}
					
					int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
					int [] pairedSitesStatAlign = statalignResult.pairedSites;
					int [] pairedSitesStatAlignWeighted = statalignWeightedResult.pairedSites;
					int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
					int [] pairedSitesMPD = mpdResult.pairedSites;
					
					
					double sensExpStat = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesStatAlign);
					double sensExpPPfold = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesPPfold);
					double sensExpMPD=Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesMPD);
					double ppvExpStat = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesStatAlign);
					double ppvExpPPfold = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesPPfold);
					double ppvExpMPD = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesMPD);
					double fscExpStat = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlign);
					double fscExpPPfold = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesPPfold);
					double fscExpStatWeighted =Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlignWeighted);
					double fscExpStatMPD=Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesMPD);
					
					System.out.println(textline+"\t"+fscExpStat+"\t"+fscExpStatMPD+"\t"+sensExpStat+"\t"+sensExpMPD+"\t"+ppvExpStat+"\t"+ppvExpMPD);
			}*/
	}
	
	public static double getDouble(double val)
	{
		if(Double.isNaN(val))
		{
			return 0;
		}
		
		return val;
	}
	
	public static void automatedTests()
	{
		String dir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/";
		String resultsDir = "/home/michael/workspace/StatAlignExecute/combined/";
		File outFile = new File("Benchmarks.txt");
		String suffix = "_5seqs";
		if(!suffix.equals("_5seqs"))
		{
			dir = "/home/michael/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		}

		
		String header = "dataset\tposterior_avg\tsim_mpd_ref\t"
		+"average_length\t"
		+"fsc_sample_mean\tfsc_sample_median\tfsc_stat\tfsc_stat_weighted\tfsc_mpd\tfsc_ppfold\tfsc_entropy_exp\tfsc_entropy_obs\t"
		+"rel_stat\trel_stat_weighted\trel_mpd\trel_ppfold\trel_entropy_exp\trel_entropy_obs\t"
		+"rel2_stat\trel2_stat_weighted\trel2_mpd\trel2_ppfold\trel2_entropy_exp\trel2_entropy_obs\t"
		+"entropy_exp\tentropy_perc_exp\tentropy_max_exp\t"
		+"entropy_obs\tentropy_perc_obs\tentropy_max_obs\t"
		+"entropy_mpd\tentropy_perc_mpd\tentropy_max_mpd\t"
		+"entropy_ppfold\tentropy_perc_ppfold\tentropy_max_ppfold\t"
		+"entropy_sample_mean\tentropy_perc_sample_mean\tentropy_max_sample_mean\t"
		+"fsc_sample_alifold_mean\tfsc_sample_alifold_median\tfsc_alifold\tfsc_alifold_mpd\talifold_ref\trel3_stat\trel3_mpd\trel3_entropy_obs";
		RNAFoldingTools.writeToFile(outFile, header, false);
		
		File [] files = new File(resultsDir).listFiles();
		for(int i = 0 ; i < files.length ; i++)
		{
			if(!files[i].getName().contains("480298957") || !files[i].getName().contains(suffix))
			{
				continue;
			}
			
			if(!files[i].getName().endsWith(".serialized"))
			{
				continue;
			}
			
			Dataset dataset = null;
			try
			{
				System.out.println(files[i]);
				dataset = Dataset.loadDatasetResult(files[i]);
			}
			catch(Exception ex)
			{
				System.err.println("Should delete "+files[i]);
				continue;
			}
			
			
			//String resultsDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq2/";
			//String dir = "/home/michael/Dropbox/RNA and StatAlign/Distance/Datasets2/";
			
			System.out.println(dataset.title);
				
			String name = dataset.title.replaceAll("_seed.+", "");
			String smallname = name.substring(0, name.length()-8);
			System.out.println(name);
			
			//String truncname = smallname;
			String truncname = smallname.replaceAll(suffix, "");
			//File files = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/").listFiles();
			File experimentalFile = new File(dir+truncname+".dat");
			//File ourData = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNADATA1OURS");
			//File statalignResultFile = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNAData1.dat.txt");
			File ppfoldData = new File(dir+truncname+".dat.ct");
			
			if(!suffix.equals("_5seqs"))
			{
				experimentalFile = new File(dir+truncname+suffix+".dat");
				ppfoldData = new File(dir+truncname+suffix+".dat.ct");
			}
			
			ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
			//StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res"));
			//StatAlignResult statalignWeightedResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res.weighted"));
			System.out.println(name);
			//StatAlignResult mpdResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res.mpd"));
			
			String mappingSeq = "";
			for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
			{
				if(experimentalData.sequences.get(j).replaceAll("-", "").equals(dataset.pairedSitesRefSeq.replaceAll("-", "")))
				{
					mappingSeq = experimentalData.sequences.get(j);
				}
			}
			
			//System.out.println(mappingSeq);
			//System.out.println(statalignResult.sequence);
			
			int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
			int [] pairedSitesStatAlign = dataset.pairedSites;
			int [] pairedSitesStatAlignWeighted = dataset.pairedSitesWeighted;
			int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
			int [] pairedSitesMPD = dataset.pairedSitesMPD;
			//printPairs(pairedSitesExperimental);
			//printPairs(pairedSitesStatAlign);
			//printPairs(pairedSitesPPfold);
			
	
			
			//System.out.println("X:"+statalignResult.sequence);
			//System.out.println("X:"+experimentalData.sequences.get(2));
			System.out.println(">"+name + " (" + (dataset.pairedSitesRefSeq.replaceAll("-", "").length())+")");
			System.out.println("  "+dataset.pairedSitesRefSeq.replaceAll("-", ""));
			System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
			System.out.println("S:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesStatAlign));
			System.out.println("  "+dataset.pairedSitesRefSeq.replaceAll("-", ""));
			System.out.println("E:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesExperimental));
			System.out.println("P:"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSitesPPfold));
			System.out.println(dataset.posteriorsAverage+"\t"+dataset.mpdVsInputSim);
			
			double sensExpStat = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesStatAlign);
			double sensExpPPfold = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesPPfold);
			double ppvExpStat = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesStatAlign);
			double ppvExpPPfold = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesPPfold);
			double fscExpStat = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlign));
			double fscExpPPfold = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesPPfold));
			double fscExpStatWeighted = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlignWeighted));
			double fscExpStatMPD= getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesMPD));
			double fscSamplingExp = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesEntropyExp));
			double fscSamplingObs = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesEntropyObs));
			double fscRNAalifold = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesRNAalifold));
			double fscCombined = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesCombined));
			System.out.println(dataset.title);
			double fscRNAalifoldMPD = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesRNAalifoldMPDProjected));
			double fscRNAalifoldRef = -1;
			if(dataset.pairedSitesRNAalifoldRefProjected != null)
			{
				fscRNAalifoldRef = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesRNAalifoldRefProjected));
			}
			else
			{
				System.err.println("Could not load.");
			}
					
			System.out.println(RNAFoldingTools.pad("", 14)+RNAFoldingTools.pad("StatAl", 10)+RNAFoldingTools.pad("PPfold", 10));
			System.out.println(RNAFoldingTools.pad("Senstivity", 14)+RNAFoldingTools.pad(sensExpStat+"", 6)+"    "+RNAFoldingTools.pad(sensExpPPfold+"", 6));
			System.out.println(RNAFoldingTools.pad("PPV", 14)+RNAFoldingTools.pad(ppvExpStat+"", 6)+"    "+RNAFoldingTools.pad(ppvExpPPfold+"", 6));
			System.out.println(RNAFoldingTools.pad("F-score", 14)+RNAFoldingTools.pad(fscExpStat+"", 6)+"    "+RNAFoldingTools.pad(fscExpPPfold+"", 6));
			System.out.print("StatAl: ");
			printValues(pairedSitesExperimental, pairedSitesStatAlign);
			System.out.print("PPfold: ");
			printValues(pairedSitesExperimental, pairedSitesPPfold);
			System.out.println("---------------------------------------------------------------------");
			System.out.println();
			ArrayList<Double> ppfoldValues = new ArrayList<Double>();
			ArrayList<Double> rnaAlifoldValues = new ArrayList<Double>();
			for(int k= 0 ; k < dataset.pairedSitesProjectedSamples.size() ; k++)
			{
				String val = "" + getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesProjectedSamples.get(k)));
				ppfoldValues.add(new Double(val));
				
			}
			
			for(int k= 0 ; k < dataset.pairedSitesProjectedRnaAlifoldSamples.size() ; k++)
			{
				//System.out.println(dataset.pairedSitesProjectedSamples.size()+"\t"+dataset.pairedSitesProjectedRnaAlifoldSamples.size());
				rnaAlifoldValues.add(getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesProjectedRnaAlifoldSamples.get(k))));
			}
			double fscSampleMean = mean(ppfoldValues);
			double fscSampleMedian = getMedian(ppfoldValues);
			
			double fscRnaAlifoldSampleMean = mean(rnaAlifoldValues);
			double fscRnaAlifoldSampleMedian = getMedian(rnaAlifoldValues);
			
			ArrayList<Double> entropySamples = new ArrayList<Double>();
			ArrayList<Double> entropyPercSamples = new ArrayList<Double>();
			ArrayList<Double> entropyMaxSamples = new ArrayList<Double>();
			for(int k = 0 ; k < dataset.sampledStructures.size() ; k++)
			{
				ResultBundle sampleBundle = dataset.sampledStructures.get(k);
				entropySamples.add(sampleBundle.entropyVal);
				entropyPercSamples.add(sampleBundle.entropyPercOfMax);
				entropyMaxSamples.add(sampleBundle.entropyMax);
			}
			double entropySampleMean = mean(entropySamples);
			double entropyPercSampleMean = mean(entropyPercSamples);
			double entropyMaxSampleMean = mean(entropyMaxSamples);


			double averageLength = getAverageLength(dataset.inputAlignment.sequences);
			String row = name+"\t"+dataset.posteriorsAverage+"\t"+dataset.mpdVsInputSim+"\t"
			+averageLength+"\t"
			+fscSampleMean+"\t"+fscSampleMedian+"\t"+fscExpStat+"\t"+fscExpStatWeighted+"\t"+fscExpStatMPD+"\t"+fscExpPPfold+"\t"+fscSamplingExp+"\t"+fscSamplingObs
			+"\t"+dataset.ppfoldReliabilityScoreSamplingAndAveraging+"\t"+dataset.ppfoldReliabilityScoreSamplingAndAveragingWeighted+"\t"+dataset.ppfoldReliabilityMPD
			+"\t"+dataset.resultBundlePPfold.reliabilityScore+"\t"+dataset.resultBundleEntropyExp.reliabilityScore+"\t"+dataset.resultBundleEntropyObs.reliabilityScore
			+"\t"+dataset.pairsOnlyReliabilityScoreSamplingAndAveraging+"\t"+dataset.pairsOnlyReliabilityScoreSamplingAndAveragingWeighted+"\t"+dataset.pairsOnlyReliabilityMPD
			+"\t"+dataset.resultBundlePPfold.pairsOnlyReliabilityScore+"\t"+dataset.resultBundleEntropyExp.pairsOnlyReliabilityScore+"\t"+dataset.resultBundleEntropyObs.pairsOnlyReliabilityScore
			+"\t"+dataset.resultBundleEntropyExp.entropyVal+"\t"+(dataset.resultBundleEntropyExp.entropyPercOfMax/100)+"\t"+dataset.resultBundleEntropyExp.entropyMax
			+"\t"+dataset.resultBundleEntropyObs.entropyVal+"\t"+(dataset.resultBundleEntropyObs.entropyPercOfMax/100)+"\t"+dataset.resultBundleEntropyObs.entropyMax
			+"\t"+dataset.resultBundleMPD.entropyVal+"\t"+(dataset.resultBundleMPD.entropyPercOfMax/100)+"\t"+dataset.resultBundleMPD.entropyMax
			+"\t"+dataset.resultBundlePPfold.entropyVal+"\t"+(dataset.resultBundlePPfold.entropyPercOfMax/100)+"\t"+dataset.resultBundlePPfold.entropyMax
			+"\t"+entropySampleMean+"\t"+(entropyPercSampleMean/100)+"\t"+entropyMaxSampleMean+"\t"
			+fscRnaAlifoldSampleMean+"\t"+fscRnaAlifoldSampleMedian+"\t"+fscRNAalifold+"\t"+fscRNAalifoldMPD+"\t"+fscRNAalifoldRef
			+"\t"+dataset.pairsOnlyReliabilityScoreSamplingAndAveragingPosteriorWeighted+"\t"+dataset.pairsOnlyMPDPosteriorWeighted+"\t"+dataset.pairsOnlyReliabilityEntropyObsPosteriorWeighted
			+"\t"+fscCombined+"\t"+dataset.pairsOnlyReliabilityScoreCombined+"\t"+dataset.ppfoldReliabilityScoreCombined;

			
		
			ArrayList<Double> fuzzyDistances = new ArrayList<Double>();
			for(int k = 1 ; k  < dataset.cumulativeFuzzyAlignment.size() ; k++)
			{
				fuzzyDistances.add(FuzzyAlignment.distance(dataset.cumulativeFuzzyAlignment.get(k-1), dataset.cumulativeFuzzyAlignment.get(k)));
			}
			//System.out.println(dataset.title+"\t"+fuzzyDistances);
			try
			{
				BufferedWriter buffer = new BufferedWriter(new FileWriter("distances/"+dataset.title+"_fuzzy_distances.txt"));
				for(int k = 0 ; k  < fuzzyDistances.size() ; k++)
				{
					buffer.write(fuzzyDistances.get(k)+"\n");	
				}
				buffer.close();				
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}
			
			
			RNAFoldingTools.writeToFile(outFile, row, true);
			
			//System.out.println(header);
			//System.out.println(row);
	
			try
			{
				System.out.println("Writing"+resultsDir+name+".hist2");
				BufferedWriter buffer = new BufferedWriter(new FileWriter(resultsDir+name+".hist2"));
				buffer.write("ST="+fscExpStat+"\n");
				buffer.write("STW="+fscExpStatWeighted+"\n");
				buffer.write("MPD="+fscExpStatMPD+"\n");
				buffer.write("PP="+fscExpPPfold+"\n");
				buffer.write("STE="+fscSamplingObs+"\n");
				for(int l = 0 ; l < ppfoldValues.size() ; l++)
				{
					double val =ppfoldValues.get(l);
					if(Double.isNaN(val))
					{
						val = 0;
					}
					buffer.write(val+"\n");
				}
				buffer.close();
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}
			
			try
			{
				File entropyFile = new File("/home/michael/Dropbox/RNA and StatAlign/Report/Entropy2/"+name+".txt");
				BufferedWriter buffer = new BufferedWriter(new FileWriter(entropyFile));
				buffer.write("no\tobs_val\tobs_perc\tobs_max\texp_val\texp_perc\texp_max\tsam_val\tsam_perc\tsam_max\n");
				for(int l = 0 ; l < dataset.sampledStructures.size() ; l++)
				{
					ResultBundle sample = dataset.sampledStructures.get(l);
					ResultBundle obsSample = dataset.cumulativeFuzzyObsResults.get(l);
					ResultBundle expSample = dataset.cumulativeFuzzyExpResults.get(l);
					buffer.write(l+"\t"+obsSample.entropyVal+"\t"+obsSample.entropyPercOfMax+"\t"+obsSample.entropyMax+"\t");
					buffer.write(expSample.entropyVal+"\t"+expSample.entropyPercOfMax+"\t"+expSample.entropyMax+"\t");
					buffer.write(sample.entropyVal+"\t"+sample.entropyPercOfMax+"\t"+sample.entropyMax);
					buffer.newLine();
				}
				buffer.close();
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}
		}
	}
	
	public static double getAverageLength(ArrayList<String> sequences)
	{
		double sum = 0;
		for(int i = 0 ; i < sequences.size() ; i++)
		{
			sum += sequences.get(i).replaceAll("-", "").length();
		}
		return sum / ((double)sequences.size());
	}
	
	public static void testVariation2()
	{
		String dir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/";
		String resultsDir = "/home/michael/workspace/StatAlignExecute/output4/";
		//File outFile = new File("Benchmarks.txt");
		String suffix = "_5seqs";
		if(!suffix.equals("_5seqs"))
		{
			dir = "/home/michael/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		}

		
		String header = "dataset\tposterior_avg\tsim_mpd_ref\t"
		+"fsc_sample_mean\tfsc_stat\tfsc_stat_weighted\tfsc_mpd\tfsc_ppfold\tfsc_entropy_exp\tfsc_entropy_obs\t"
		+"rel_stat\trel_stat_weighted\trel_mpd\trel_ppfold\trel_entropy_exp\trel_entropy_obs\t"
		+"rel2_stat\trel2_stat_weighted\trel2_mpd\trel2_ppfold\trel2_entropy_exp\trel2_entropy_obs\t"
		+"entropy_exp\tentropy_perc_exp\tentropy_max_exp\t"
		+"entropy_obs\tentropy_perc_obs\tentropy_max_obs\t"
		+"entropy_mpd\tentropy_perc_mpd\tentropy_max_mpd\t"
		+"entropy_ppfold\tentropy_perc_ppfold\tentropy_max_ppfold\t"
		+"entropy_sample_mean\tentropy_perc_sample_mean\tentropy_max_sample_mean\t";
		//RNAFoldingTools.writeToFile(outFile, header, false);
		
		File [] files = new File(resultsDir).listFiles();
		HashSet<String> usedDatasets = new HashSet<String>();
		for(int z = 0 ; z < files.length ; z++)
		{
			String datasetName = files[z].getName().replaceAll("_seed.*", "");
			if(!files[z].getName().contains(suffix) || !files[z].getName().endsWith(".serialized") || usedDatasets.contains(datasetName))
			{
				continue;
			}
			
			usedDatasets.add(datasetName);
			
			ArrayList<Double> samplesFscVector = new ArrayList<Double>();
			ArrayList<Double> statalignFscVector = new ArrayList<Double>();
			ArrayList<Double> mpdFscVector = new ArrayList<Double>();
			ArrayList<Double> ppfoldFscVector = new ArrayList<Double>();
			ArrayList<Double> entropyObsFscVector = new ArrayList<Double>();
			ArrayList<Double> entropyExpFscVector = new ArrayList<Double>();
			for(int i = 0 ; i < files.length ; i++)
			{				
				if(!files[i].getName().contains(datasetName) || !files[i].getName().contains(suffix) || !files[i].getName().endsWith(".serialized"))
				{
					continue;
				}
				//System.out.println(files[i]);
				
				if(files[i].getName().equals("TestRNAData27_5seqs.dat.fas_seed682981838.serialized"))
				{
					continue;
				}
				
				System.out.println(files[i]);
				Dataset dataset = Dataset.loadDatasetResult(files[i]);
				
				
				
				
				//String resultsDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq2/";
				//String dir = "/home/michael/Dropbox/RNA and StatAlign/Distance/Datasets2/";
				
				//System.out.println(dataset.title);
					
				String name = dataset.title.replaceAll("_seed.+", "");
				String smallname = name.substring(0, name.length()-8);
				//System.out.println(name);
				
				//String truncname = smallname;
				String truncname = smallname.replaceAll(suffix, "");
				//File files = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/").listFiles();
				File experimentalFile = new File(dir+truncname+".dat");
				//File ourData = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNADATA1OURS");
				//File statalignResultFile = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNAData1.dat.txt");
				File ppfoldData = new File(dir+truncname+".dat.ct");
				
				if(!suffix.equals("_5seqs"))
				{
					experimentalFile = new File(dir+truncname+suffix+".dat");
					ppfoldData = new File(dir+truncname+suffix+".dat.ct");
				}
				
				ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
				//StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res"));
				//StatAlignResult statalignWeightedResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res.weighted"));
				//System.out.println(name);
				//StatAlignResult mpdResult = loadStatAlignResultFile(new File(resultsDir+"/"+name+".dat.res.mpd"));
				
				String mappingSeq = "";
				for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
				{
					if(experimentalData.sequences.get(j).replaceAll("-", "").equals(dataset.pairedSitesRefSeq.replaceAll("-", "")))
					{
						mappingSeq = experimentalData.sequences.get(j);
					}
				}
				
				//System.out.println(mappingSeq);
				//System.out.println(statalignResult.sequence);
				
				int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
				int [] pairedSitesStatAlign = dataset.pairedSites;
				int [] pairedSitesStatAlignWeighted = dataset.pairedSitesWeighted;
				int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
				int [] pairedSitesMPD = dataset.pairedSitesMPD;
					
				double sensExpStat = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesStatAlign);
				double sensExpPPfold = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesPPfold);
				double ppvExpStat = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesStatAlign);
				double ppvExpPPfold = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesPPfold);
				double fscExpStat = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlign));
				double fscExpPPfold = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesPPfold));
				double fscExpStatWeighted = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlignWeighted));
				double fscExpStatMPD= getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesMPD));
				double fscSamplingExp = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesEntropyExp));
				double fscSamplingObs = getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesEntropyObs));
						
				ArrayList<Double> values = new ArrayList<Double>();
				for(int k= 0 ; k < dataset.pairedSitesProjectedSamples.size() ; k++)
				{
					String val = "" + getDouble(Benchmarks.calculateFScore(pairedSitesExperimental, dataset.pairedSitesProjectedSamples.get(k)));
					values.add(new Double(val));
				}
				double fscSampleMean = mean(values);
				
				samplesFscVector.addAll(values);
				statalignFscVector.add(fscExpStat);
				mpdFscVector.add(fscExpStatMPD);
				ppfoldFscVector.add(fscExpPPfold);
				entropyObsFscVector.add(fscSamplingObs);
				entropyExpFscVector.add(fscSamplingExp);
				
				ArrayList<Double> entropySamples = new ArrayList<Double>();
				ArrayList<Double> entropyPercSamples = new ArrayList<Double>();
				ArrayList<Double> entropyMaxSamples = new ArrayList<Double>();
				for(int k = 0 ; k < dataset.sampledStructures.size() ; k++)
				{
					ResultBundle sampleBundle = dataset.sampledStructures.get(k);
					entropySamples.add(sampleBundle.entropyVal);
					entropyPercSamples.add(sampleBundle.entropyPercOfMax);
					entropyMaxSamples.add(sampleBundle.entropyMax);
				}
				double entropySampleMean = mean(entropySamples);
				double entropyPercSampleMean = mean(entropyPercSamples);
				double entropyMaxSampleMean = mean(entropyMaxSamples);
				
				String row = name+"\t"+dataset.posteriorsAverage+"\t"+dataset.mpdVsInputSim+"\t"
				+fscSampleMean+"\t"+fscExpStat+"\t"+fscExpStatWeighted+"\t"+fscExpStatMPD+"\t"+fscExpPPfold+"\t"+fscSamplingExp+"\t"+fscSamplingObs
				+"\t"+dataset.ppfoldReliabilityScoreSamplingAndAveraging+"\t"+dataset.ppfoldReliabilityScoreSamplingAndAveragingWeighted+"\t"+dataset.ppfoldReliabilityMPD
				+"\t"+dataset.resultBundlePPfold.reliabilityScore+"\t"+dataset.resultBundleEntropyExp.reliabilityScore+"\t"+dataset.resultBundleEntropyObs.reliabilityScore
				+"\t"+dataset.pairsOnlyReliabilityScoreSamplingAndAveraging+"\t"+dataset.pairsOnlyReliabilityScoreSamplingAndAveragingWeighted+"\t"+dataset.pairsOnlyReliabilityMPD
				+"\t"+dataset.resultBundlePPfold.pairsOnlyReliabilityScore+"\t"+dataset.resultBundleEntropyExp.pairsOnlyReliabilityScore+"\t"+dataset.resultBundleEntropyObs.pairsOnlyReliabilityScore
				+"\t"+dataset.resultBundleEntropyExp.entropyVal+"\t"+(dataset.resultBundleEntropyExp.entropyPercOfMax/100)+"\t"+dataset.resultBundleEntropyExp.entropyMax
				+"\t"+dataset.resultBundleEntropyObs.entropyVal+"\t"+(dataset.resultBundleEntropyObs.entropyPercOfMax/100)+"\t"+dataset.resultBundleEntropyObs.entropyMax
				+"\t"+dataset.resultBundleMPD.entropyVal+"\t"+(dataset.resultBundleMPD.entropyPercOfMax/100)+"\t"+dataset.resultBundleMPD.entropyMax
				+"\t"+dataset.resultBundlePPfold.entropyVal+"\t"+(dataset.resultBundlePPfold.entropyPercOfMax/100)+"\t"+dataset.resultBundlePPfold.entropyMax
				+"\t"+entropySampleMean+"\t"+(entropyPercSampleMean/100)+"\t"+entropyMaxSampleMean;
			
				/*
				ArrayList<Double> fuzzyDistances = new ArrayList<Double>();
				for(int k = 1 ; k  < dataset.cumulativeFuzzyAlignment.size() ; k++)
				{
					fuzzyDistances.add(FuzzyAlignment.distance(dataset.cumulativeFuzzyAlignment.get(k-1), dataset.cumulativeFuzzyAlignment.get(k)));
				}
				//System.out.println(dataset.title+"\t"+fuzzyDistances);
				try
				{
					BufferedWriter buffer = new BufferedWriter(new FileWriter("distances/"+dataset.title+"_fuzzy_distances.txt"));
					for(int k = 0 ; k  < fuzzyDistances.size() ; k++)
					{
						buffer.write(fuzzyDistances.get(k)+"\n");	
					}
					buffer.close();				
				}
				catch(IOException ex)
				{
					ex.printStackTrace();
				}*/
			}
			
			DecimalFormat df = new DecimalFormat("0.000");
			//System.out.println();
			String vdir = "/home/michael/Dropbox/RNA and StatAlign/Report/V2/";
			//System.out.println("Dataset\t\t\t\t\t#\tsample\tstat\tmpd\tppfold\ten_exp\ten_obs");
			RNAFoldingTools.writeToFile(new File(vdir+"mean.txt"), datasetName+"\tmean\t\t"+mpdFscVector.size()+"\t"+df.format(mean(samplesFscVector))+"\t"+df.format(mean(statalignFscVector))+"\t"+df.format(mean(mpdFscVector))+"\t"+df.format(mean(ppfoldFscVector))+"\t"+df.format(mean(entropyExpFscVector))+"\t"+df.format(mean(entropyObsFscVector)), true);
			//System.out.println(datasetName+"\tmean\t\t"+mpdFscVector.size()+"\t"+df.format(mean(samplesFscVector))+"\t"+df.format(mean(statalignFscVector))+"\t"+df.format(mean(mpdFscVector))+"\t"+df.format(mean(ppfoldFscVector))+"\t"+df.format(mean(entropyExpFscVector))+"\t"+df.format(mean(entropyObsFscVector)));
			RNAFoldingTools.writeToFile(new File(vdir+"stdev.txt"), datasetName+"\tstdev\t\t"+mpdFscVector.size()+"\t"+df.format(stdev(samplesFscVector))+"\t"+df.format(stdev(statalignFscVector))+"\t"+df.format(stdev(mpdFscVector))+"\t"+df.format(stdev(ppfoldFscVector))+"\t"+df.format(stdev(entropyExpFscVector))+"\t"+df.format(stdev(entropyObsFscVector)), true);
			//System.out.println(datasetName+"\tstdev\t\t"+mpdFscVector.size()+"\t"+df.format(stdev(samplesFscVector))+"\t"+df.format(stdev(statalignFscVector))+"\t"+df.format(stdev(mpdFscVector))+"\t"+df.format(stdev(ppfoldFscVector))+"\t"+df.format(stdev(entropyExpFscVector))+"\t"+df.format(stdev(entropyObsFscVector)));
			double perc = 0.25;
			RNAFoldingTools.writeToFile(new File(vdir+"25th.txt"), datasetName+"\t25th\t\t"+mpdFscVector.size()+"\t"+df.format(getValue(samplesFscVector, perc))+"\t"+df.format(getValue(statalignFscVector, perc))+"\t"+df.format(getValue(mpdFscVector, perc))+"\t"+df.format(getValue(ppfoldFscVector, perc))+"\t"+df.format(getValue(entropyExpFscVector, perc))+"\t"+df.format(getValue(entropyObsFscVector, perc)), true);
			//System.out.println(datasetName+"\t25th\t\t"+mpdFscVector.size()+"\t"+df.format(getValue(samplesFscVector, perc))+"\t"+df.format(getValue(statalignFscVector, perc))+"\t"+df.format(getValue(mpdFscVector, perc))+"\t"+df.format(getValue(ppfoldFscVector, perc))+"\t"+df.format(getValue(entropyExpFscVector, perc))+"\t"+df.format(getValue(entropyObsFscVector, perc)));
			perc = 0.5;
			RNAFoldingTools.writeToFile(new File(vdir+"50th.txt"), datasetName+"\tmedian\t\t"+mpdFscVector.size()+"\t"+df.format(getValue(samplesFscVector, perc))+"\t"+df.format(getValue(statalignFscVector, perc))+"\t"+df.format(getValue(mpdFscVector, perc))+"\t"+df.format(getValue(ppfoldFscVector, perc))+"\t"+df.format(getValue(entropyExpFscVector, perc))+"\t"+df.format(getValue(entropyObsFscVector, perc)), true);
			//System.out.println(datasetName+"\tmedian\t\t"+mpdFscVector.size()+"\t"+df.format(getValue(samplesFscVector, perc))+"\t"+df.format(getValue(statalignFscVector, perc))+"\t"+df.format(getValue(mpdFscVector, perc))+"\t"+df.format(getValue(ppfoldFscVector, perc))+"\t"+df.format(getValue(entropyExpFscVector, perc))+"\t"+df.format(getValue(entropyObsFscVector, perc)));
			perc = 0.75;
			RNAFoldingTools.writeToFile(new File(vdir+"75th.txt"),datasetName+"\t75th\t\t"+mpdFscVector.size()+"\t"+df.format(getValue(samplesFscVector, perc))+"\t"+df.format(getValue(statalignFscVector, perc))+"\t"+df.format(getValue(mpdFscVector, perc))+"\t"+df.format(getValue(ppfoldFscVector, perc))+"\t"+df.format(getValue(entropyExpFscVector, perc))+"\t"+df.format(getValue(entropyObsFscVector, perc)), true);
			//System.out.println(datasetName+"\t75th\t\t"+mpdFscVector.size()+"\t"+df.format(getValue(samplesFscVector, perc))+"\t"+df.format(getValue(statalignFscVector, perc))+"\t"+df.format(getValue(mpdFscVector, perc))+"\t"+df.format(getValue(ppfoldFscVector, perc))+"\t"+df.format(getValue(entropyExpFscVector, perc))+"\t"+df.format(getValue(entropyObsFscVector, perc)));
			RNAFoldingTools.writeToFile(new File(vdir+"IQR.txt"),datasetName+"\tIQR\t\t"+mpdFscVector.size()+"\t"+df.format(IQR(samplesFscVector))+"\t"+df.format(IQR(statalignFscVector))+"\t"+df.format(IQR(mpdFscVector))+"\t"+df.format(IQR(ppfoldFscVector))+"\t"+df.format(IQR(entropyExpFscVector))+"\t"+df.format(IQR(entropyObsFscVector)), true);
			System.out.println(datasetName+"\tIQR\t\t"+mpdFscVector.size()+"\t"+df.format(IQR(samplesFscVector))+"\t"+df.format(IQR(statalignFscVector))+"\t"+df.format(IQR(mpdFscVector))+"\t"+df.format(IQR(ppfoldFscVector))+"\t"+df.format(IQR(entropyExpFscVector))+"\t"+df.format(IQR(entropyObsFscVector)));
			double sampleMean = mean(samplesFscVector);
			RNAFoldingTools.writeToFile(new File(vdir+"percent_greater_than_mean.txt"),datasetName+"\t%>mean\t\t"+mpdFscVector.size()+"\t"+df.format(percentGreaterOrEqualTo(sampleMean, samplesFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean, statalignFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean,mpdFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean,ppfoldFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean,entropyExpFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean,entropyObsFscVector)), true);
			double sampleMedian = getValue(samplesFscVector, 0.5);
			RNAFoldingTools.writeToFile(new File(vdir+"percent_greater_than_median.txt"),datasetName+"\t%>median\t\t"+mpdFscVector.size()+"\t"+df.format(percentGreaterOrEqualTo(sampleMedian, samplesFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMedian, statalignFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMedian,mpdFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMedian,ppfoldFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMedian,entropyExpFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMedian,entropyObsFscVector)), true);
			
			try
			{
				BufferedWriter buffer = new BufferedWriter(new FileWriter(System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/Report/Variation/"+datasetName+".var"));
				buffer.write("Samples\tStatAlign\tMPD\tPPfold\tEntropy obs\n");
				for(int k = 0 ; k < statalignFscVector.size() ; k++)
				{
					buffer.write(samplesFscVector.get(k)+"\t"+statalignFscVector.get(k)+"\t"+mpdFscVector.get(k)+"\t"+ppfoldFscVector.get(k)+"\t"+entropyObsFscVector.get(k));
					buffer.newLine();
				}
				for(int k = statalignFscVector.size() ; k < samplesFscVector.size() ; k++)
				{
					buffer.write(samplesFscVector.get(k)+"\n");
				}
				buffer.close();
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}
		}
	}
	
	public static void testVariation()
	{
		
		String dir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		//String dir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/";
		//String resultsDir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/Results3/";
		String resultsDir = "/home/michael/workspace/StatAlignExecute/output/";
		
		File [] files = new File(resultsDir).listFiles();
		HashSet<String> used = new HashSet<String>();
		for(int i = 0 ; i < files.length ; i++)
		{
			String fullName = files[i].getName();
			if(fullName.endsWith(".dat.res"))
			{
				String name = fullName.substring(0, fullName.length()-8);
				String smallname = fullName.substring(0, fullName.length()-16);
				//System.out.println(fullName);
				//System.out.println(name);
				//System.out.println(smallname);
		
				String originalName = name.replaceAll("_seed.+", "");
				originalName = originalName.substring(0, originalName.length()-8);
				//System.out.println("O:"+originalName);
				File experimentalFile = new File(dir+originalName+".dat");			
				File ppfoldData = new File(dir+originalName+".dat.ct");
				
				ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
				

				if(used.contains(originalName))
				{
					continue;
				}
				used.add(originalName);
							
				int runs = 0;
				ArrayList<Double> samplesFscVector = new ArrayList<Double>();
				ArrayList<Double> statalignFscVector = new ArrayList<Double>();
				ArrayList<Double> mpdFscVector = new ArrayList<Double>();
				ArrayList<Double> ppfoldFscVector = new ArrayList<Double>();
				ArrayList<Double> entropyObsFscVector = new ArrayList<Double>();
				ArrayList<Double> entropyExpFscVector = new ArrayList<Double>();
				for(int l = 0; l < files.length ; l++)
				{
					
					if(files[l].getName().startsWith(originalName) && files[l].getName().endsWith(".dat.res"))
					{
						String runName = files[l].getName().substring(0, files[l].getName().length()-8);
						runs++;
					
						StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+runName+".dat.res"));
						StatAlignResult statalignWeightedResult = loadStatAlignResultFile(new File(resultsDir+"/"+runName+".dat.res.weighted"));
						if(!new File(resultsDir+"/"+runName+".dat.res.mpd").exists())
						{
							//System.err.println(new File(resultsDir+"/"+runName+".dat.res.mpd")+" is missing");
							continue;
						}
						StatAlignResult mpdResult = loadStatAlignResultFile(new File(resultsDir+"/"+runName+".dat.res.mpd"));
						
						
						String mappingSeq = "";
						for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
						{
							if(experimentalData.sequences.get(j).replaceAll("-", "").equals(statalignResult.sequence.replaceAll("-", "")))
							{
								mappingSeq = experimentalData.sequences.get(j);
							}
						}
						
						int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
						int [] pairedSitesStatAlign = statalignResult.pairedSites;
						int [] pairedSitesStatAlignWeighted = statalignWeightedResult.pairedSites;
						int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
						int [] pairedSitesMPD = mpdResult.pairedSites;
								
					
						
						Scores statalignScores = Scores.getScores(pairedSitesExperimental, pairedSitesStatAlign);
						Scores statalignWeightedScores = Scores.getScores(pairedSitesExperimental, pairedSitesStatAlignWeighted);
						Scores ppfoldScores = Scores.getScores(pairedSitesExperimental, pairedSitesPPfold);
						Scores mpdScores = Scores.getScores(pairedSitesExperimental, pairedSitesMPD);
						double fscSamplingObs = -1;
						double fscSamplingExp = -1;
						
				
						File obsFile = new File(resultsDir+runName+".folds_e_obs");
						if(obsFile.exists())
						{
							String dbn = PPFold.loadFolds(obsFile, 4).get(0);
							fscSamplingObs = Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(dbn));
							if(Double.isNaN(fscSamplingObs))
							{
								fscSamplingObs = 0;
							}
							//System.out.println("Alignment sampling obs FSC"+fscSamplingObs);
							
							dbn = PPFold.loadFolds(new File(resultsDir+runName+".folds_e_exp"), 4).get(0);
							fscSamplingExp = Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(dbn));
							if(Double.isNaN(fscSamplingExp))
							{
								fscSamplingExp = 0;
							}
							//System.out.println("Alignment sampling exp FSC"+fscSamplingExp);
							
						}
		
		
						if(new File(resultsDir+runName+".folds").exists() && fscSamplingObs != -1)
						{
							ArrayList<String> structures = PPFold.loadFolds(new File(resultsDir+runName+".folds"), 4);
							ArrayList<Double> sampleValues = new ArrayList<Double>();
							for(int k= 0 ; k < structures.size() ; k++)
							{
								String val = "" + Benchmarks.calculateFScore(pairedSitesExperimental, RNAFoldingTools.getPairedSitesFromDotBracketString(structures.get(k)));
								
								Double d = new Double(val);
								if(d.isNaN())
								{
									d = new Double(0);
								}
								sampleValues.add(d);
							}
							

							samplesFscVector.addAll(sampleValues);
							statalignFscVector.add(statalignScores.fsc);
							ppfoldFscVector.add(ppfoldScores.fsc);
							mpdFscVector.add(mpdScores.fsc);
							entropyObsFscVector.add(fscSamplingObs);
							entropyExpFscVector.add(fscSamplingExp);
							//System.out.println(runName+"\t"+runs+"\t"+mean(values)+"\t"+statalignScores.fsc+"\t"+mpdScores.fsc+"\t"+ppfoldScores.fsc+"\t"+fscSamplingObs);
		
							/*try
							{
								BufferedWriter buffer = new BufferedWriter(new FileWriter(resultsDir+name+".hist2"));
								buffer.write("ST="+statalignScores.fsc+"\n");
								buffer.write("STW="+statalignWeightedScores.fsc+"\n");
								buffer.write("MPD="+mpdScores.fsc+"\n");
								buffer.write("PP="+ppfoldScores.fsc+"\n");
								buffer.write("STE="+fscSamplingObs+"\n");
								for(int l = 0 ; l < values.size() ; l++)
								{
									double val = Double.parseDouble(values.get(l));
									if(Double.isNaN(val))
									{
										val = 0;
									}
									buffer.write(val+"\n");
								}
								buffer.close();
							}
							catch(IOException ex)
							{
								ex.printStackTrace();
							}*/
						}
					}
				}
				
				DecimalFormat df = new DecimalFormat("0.000");
				//System.out.println();
				System.out.println("Dataset\t\t\t\t\t#\tsample\tstat\tmpd\tppfold\ten_exp\ten_obs");
				System.out.println(originalName+"\tmean\t\t"+runs+"\t"+df.format(mean(samplesFscVector))+"\t"+df.format(mean(statalignFscVector))+"\t"+df.format(mean(mpdFscVector))+"\t"+df.format(mean(ppfoldFscVector))+"\t"+df.format(mean(entropyExpFscVector))+"\t"+df.format(mean(entropyObsFscVector)));
				System.out.println(originalName+"\tstdev\t\t"+runs+"\t"+df.format(stdev(samplesFscVector))+"\t"+df.format(stdev(statalignFscVector))+"\t"+df.format(stdev(mpdFscVector))+"\t"+df.format(stdev(ppfoldFscVector))+"\t"+df.format(stdev(entropyExpFscVector))+"\t"+df.format(stdev(entropyObsFscVector)));
				double perc = 0.25;
				System.out.println(originalName+"\tstdev\t\t"+runs+"\t"+df.format(stdev(samplesFscVector))+"\t"+df.format(stdev(statalignFscVector))+"\t"+df.format(stdev(mpdFscVector))+"\t"+df.format(stdev(ppfoldFscVector))+"\t"+df.format(stdev(entropyExpFscVector))+"\t"+df.format(stdev(entropyObsFscVector)));
				double sampleMean = mean(samplesFscVector);
				
				try
				{
					BufferedWriter buffer = new BufferedWriter(new FileWriter(System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/Report/Variation/"+originalName+".var"));
					buffer.write("Samples\tStatAlign\tMPD\tPPfold\tEntropy exp\tEntropy obs\n");
					for(int k = 0 ; k < statalignFscVector.size() ; k++)
					{
						buffer.write(samplesFscVector.get(k)+"\t"+statalignFscVector.get(k)+"\t"+mpdFscVector.get(k)+"\t"+ppfoldFscVector.get(k)+"\t"+entropyExpFscVector.get(k)+"\t"+entropyObsFscVector.get(k));
						buffer.newLine();
					}
					for(int k = statalignFscVector.size() ; k < samplesFscVector.size() ; k++)
					{
						buffer.write(samplesFscVector.get(k)+"\n");
					}
					buffer.close();
				}
				catch(IOException ex)
				{
					ex.printStackTrace();
				}
				//System.out.println(originalName+"\t% >= mean\t"+runs+"\t"+df.format(percentGreaterOrEqualTo(sampleMean, samplesFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean, statalignFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean, mpdFscVector))+"\t"+df.format(score(sampleMean, ppfoldFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean, entropyExpFscVector))+"\t"+df.format(percentGreaterOrEqualTo(sampleMean, entropyObsFscVector)));
				//System.out.println(originalName+"\tscore\t\t"+runs+"\t"+df.format(score(sampleMean, samplesFscVector))+"\t"+df.format(score(sampleMean, statalignFscVector))+"\t"+df.format(score(sampleMean, mpdFscVector))+"\t"+df.format(score(sampleMean, ppfoldFscVector))+"\t"+df.format(score(sampleMean, entropyExpFscVector))+"\t"+df.format(score(sampleMean, entropyObsFscVector)));
			}
		}
	}
	
	public void performDistanceBenchmarks()
	{
		File distanceFile = new File(System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/9seq2/dist_scores.txt");
		
		String dataDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		String resultsDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/9seq2/";
		
		//File distanceFile = new File(System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq/dist_scores.txt");
		
		//String dataDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/TestRNAData/";
		//String resultsDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq/";
		try
		{
			BufferedReader buffer = new BufferedReader(new FileReader(distanceFile));
			String textline = null;
			String header = "";
			while((textline = buffer.readLine()) != null)
			{
				if(textline.startsWith("Dataset"))
				{
					header = textline;
					//System.out.println(header);
					continue;
				}
				String [] split = textline.split("\t+");
				String dataName = split[0];
				String name = dataName.split("\\.")[0];
					
					
					File experimentalFile = new File(dataDir+name+".dat");
					File ppfoldData = new File(dataDir+name+".dat.ct");
					
					ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
					
					StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res")); // 
					StatAlignResult statalignWeightedResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res.weighted"));
					StatAlignResult mpdResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res.mpd"));
					
					String mappingSeq = "";
					for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
					{
						if(experimentalData.sequences.get(j).replaceAll("-", "").equals(statalignResult.sequence.replaceAll("-", "")))
						{
							mappingSeq = experimentalData.sequences.get(j);
						}
					}
					
					int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
					int [] pairedSitesStatAlign = statalignResult.pairedSites;
					int [] pairedSitesStatAlignWeighted = statalignWeightedResult.pairedSites;
					int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
					int [] pairedSitesMPD = mpdResult.pairedSites;
					
					
					double sensExpStat = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesStatAlign);
					double sensExpPPfold = Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesPPfold);
					double sensExpMPD=Benchmarks.calculateSensitivity(pairedSitesExperimental, pairedSitesMPD);
					double ppvExpStat = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesStatAlign);
					double ppvExpPPfold = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesPPfold);
					double ppvExpMPD = Benchmarks.calculatePPV(pairedSitesExperimental, pairedSitesMPD);
					double fscExpStat = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlign);
					double fscExpPPfold = Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesPPfold);
					double fscExpStatWeighted =Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlignWeighted);
					double fscExpStatMPD=Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesMPD);
					
					System.out.println(textline+"\t"+fscExpStat+"\t"+fscExpStatMPD+"\t"+sensExpStat+"\t"+sensExpMPD+"\t"+ppvExpStat+"\t"+ppvExpMPD);
			}
			
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	public void performEntropy()
	{
		File distanceFile = new File(System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq2/dist_scores.txt");
		
		//String dataDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		String resultsDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq2/";
		
		//File distanceFile = new File(System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq/dist_scores.txt");
		
		String dataDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/TestRNAData/";
		//String resultsDir = System.getProperty("user.home")+ "/Dropbox/RNA and StatAlign/static/5seq/";
		try
		{
			BufferedReader buffer = new BufferedReader(new FileReader(distanceFile));
			String textline = null;
			String header = "";
			while((textline = buffer.readLine()) != null)
			{
				if(textline.startsWith("Dataset"))
				{
					header = textline;
					//System.out.println(header);
					continue;
				}
				String [] split = textline.split("\t+");
				String dataName = split[0];
				String name = dataName.split("\\.")[0].split("_")[0];
					
					
				File experimentalFile = new File(dataDir+name+".dat");
				File ppfoldData = new File(dataDir+name+".dat.ct");
				
				ExperimentalData experimentalData = Benchmarks.loadExperimentalStructure(experimentalFile);
				
				StatAlignResult statalignResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res")); // 
				StatAlignResult statalignWeightedResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res.weighted"));
				StatAlignResult mpdResult = loadStatAlignResultFile(new File(resultsDir+"/"+dataName+".dat.res.mpd"));
				
				String mappingSeq = "";
				for(int j = 0 ;j < experimentalData.sequences.size() ; j++)
				{
					if(experimentalData.sequences.get(j).replaceAll("-", "").equals(statalignResult.sequence.replaceAll("-", "")))
					{
						mappingSeq = experimentalData.sequences.get(j);
					}
				}
				
				int [] pairedSitesExperimental = projectPairedSites(mappingSeq, experimentalData.pairedSites);
				int [] pairedSitesStatAlign = statalignResult.pairedSites;
				int [] pairedSitesStatAlignWeighted = statalignWeightedResult.pairedSites;
				int [] pairedSitesPPfold = projectPairedSites(mappingSeq, RNAFoldingTools.getPairedSitesFromCtFile(ppfoldData));
				int [] pairedSitesMPD = mpdResult.pairedSites;
				
				
				Scores statalignScores = Scores.getScores(pairedSitesExperimental, pairedSitesStatAlign);
				Scores ppfoldScores = Scores.getScores(pairedSitesExperimental, pairedSitesPPfold);
				Scores mpdScores = Scores.getScores(pairedSitesExperimental, pairedSitesMPD);
				
		
				double fscExpStatWeighted =Benchmarks.calculateFScore(pairedSitesExperimental, pairedSitesStatAlignWeighted);

				
				EntropyData entropyDataExp = EntropyData.loadEntropyData(new File(resultsDir+dataName+"_entropy_fuzzy_exp.txt"));
				EntropyData entropyDataObs = EntropyData.loadEntropyData(new File(resultsDir+dataName+"_entropy_fuzzy_obs.txt"));
				EntropyData entropyDataSamples = EntropyData.loadEntropyData(new File(resultsDir+dataName+"_entropy_samples.txt"));
				
				
				double entropyDataExpLast = entropyDataExp.entropyVals.get(entropyDataExp.entropyVals.size()-1);
				double entropyDataObsLast = entropyDataObs.entropyVals.get(entropyDataExp.entropyVals.size()-1);
				double entropyDataObsPercLast = entropyDataObs.percentOfMax.get(entropyDataExp.percentOfMax.size()-1);
				double sampleMean = mean(entropyDataSamples.entropyVals);
				double sampleMeanPercentOfMax = mean(entropyDataSamples.percentOfMax);
				
				System.out.println(textline+"\t"+statalignScores.fsc+"\t"+mpdScores.fsc+"\t"+statalignScores.sen+"\t"+mpdScores.sen+"\t"+statalignScores.ppv+"\t"+mpdScores.ppv+"\t"+entropyDataObsLast+"\t"+sampleMean+"\t"+entropyDataObsPercLast+"\t"+sampleMeanPercentOfMax);
			}
			
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	static class EntropyData
	{
		ArrayList<Double> sampleNo = new ArrayList<Double>();
		ArrayList<Double> entropyVals = new ArrayList<Double>();
		ArrayList<Double> percentOfMax = new ArrayList<Double>();
		ArrayList<Double> maxVals = new ArrayList<Double>();
		
		public static EntropyData loadEntropyData(File entropyFile)
		{
			EntropyData entropyData = new EntropyData();
			
			try
			{
				BufferedReader buffer = new BufferedReader(new FileReader(entropyFile));
				buffer.readLine();
				String textline = null;
				while((textline = buffer.readLine()) != null)
				{
					String [] split = textline.split("(\t)+");
					entropyData.sampleNo.add(new Double(split[0]));
					entropyData.entropyVals.add(new Double(split[1]));
					entropyData.percentOfMax.add(new Double(split[2]));
					entropyData.maxVals.add(new Double(split[3]));
				}
				buffer.close();
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}
			
			return entropyData;
		}
	}
	
	public static double percentGreaterOrEqualTo(double x, ArrayList<Double> values)
	{
		double count = 0;
		for(int i = 0 ; i < values.size() ; i++)
		{
			if(values.get(i) >= x)
			{
				count++;
			}
		}
		
		return count / ((double)values.size());
	}
	
	public static double score(double x, ArrayList<Double> values)
	{
		double count = 0;
		for(int i = 0 ; i < values.size() ; i++)
		{
			count += Math.min(values.get(i)/x, 1);
		}
		
		return count / ((double)values.size());
	}
	
	public static double mean(ArrayList<Double> values)
	{
		double sum = 0;
		for(int i = 0 ; i < values.size() ; i++)
		{
			sum += values.get(i);
		}
		return sum / ((double) values.size());
	}
	
	public static double stdev(ArrayList<Double> values)
	{
		double mean = mean(values);
		double stdev = 0;
		for(int i = 0 ; i < values.size() ; i++)
		{
			stdev += Math.pow(values.get(i) - mean, 2);
		}
		stdev /= ((double)values.size()-1);
		return Math.sqrt(stdev);
	}
	
	public static void printPairs(int [] pairedSites)
	{
		String ret = "";
		for(int i = 0 ; i < pairedSites.length ; i++)
		{
			ret += (i+1) + "\t" + pairedSites[i]+"\n";
		}
		System.out.println(ret);
	}
	
	public static int [] projectPairedSites(String alignedSequence, int [] pairedSites)
	{
		int [] ungappedToGapped = Mapping.getUngappedToGappedMapping(alignedSequence);
		int [] gappedToUngapped = Mapping.getGappedToUngappedMapping(alignedSequence);
		
		int [] projectedPairedSites = new int[ungappedToGapped.length];
		for(int i = 0 ; i < pairedSites.length ; i++)
		{
			if(pairedSites[i] != 0) // if paired, map
			{
				int x = gappedToUngapped[i];
				if(x != -1);
				{
					int y = Math.max(0, gappedToUngapped[pairedSites[i]-1]) + 1;
					projectedPairedSites[x] = y;
				}
			}
		}
		return projectedPairedSites;
	}
	
	public static void saveAsFasta(ExperimentalData expData, File outFile)
	{
		saveAsFasta(expData.sequences, expData.sequenceNames, outFile);
	}
	
	public static void saveAsFasta(ArrayList<String> sequences, ArrayList<String> sequenceNames, File outFile)
	{
		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile));
			for(int i = 0 ; i < sequences.size() ; i++)
			{
				buffer.write(">"+sequenceNames.get(i)+"\n");
				buffer.write(sequences.get(i)+"\n");
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	/**
	 * Loads the experimentally-derived secondary structures and alignments.
	 * @param realStructureFile the name of the file to load.
	 * @return an ExperimentalData object representing the experimental structure and alignment.
	 */
	public static ExperimentalData loadExperimentalStructure(File realStructureFile)
	{
		try
		{
			ExperimentalData expData = new ExperimentalData();
			BufferedReader buffer = new BufferedReader(new FileReader(realStructureFile));
			String realStructureDBS = buffer.readLine();
			
			expData.pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketString(realStructureDBS, '<', '>');
			ArrayList<String> sequences = new ArrayList<String>();
			ArrayList<String> sequenceNames = new ArrayList<String>();
			
			String textline = null;
            String sequence = "";
            while ((textline = buffer.readLine()) != null) {                
                if (textline.startsWith(">")) {
                    sequenceNames.add(textline.substring(1));
                    if (!sequence.equals("")) {
                        sequences.add(sequence.toUpperCase());
                        sequence = "";
                    }
                } else {
                    sequence += textline.trim();
                }

            }
            buffer.close();
            if (!sequence.equals("")) {
                sequences.add(sequence);
            }

            // replace .'s with -'s.
            for(int i = 0 ; i < sequences.size() ; i++)
            {
            	sequences.set(i, sequences.get(i).replaceAll("\\.", "-"));
            }
            
            expData.sequences = sequences;
            expData.sequenceNames = sequenceNames;
            
            return expData;			
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
		
		return null;
	}
	
	/**
	 * Given an array of paired sites corresponding to the real structure 
	 * and an array corrsponding to the predicted structure returns the sensitivity.
	 * @param realPairedSites
	 * @param predictedPairedSites
	 */
	public static double calculateSensitivity (int [] realPairedSites, int [] predictedPairedSites)
	{/* The sensitivity for a predicted structure
		is the percentage of base pairs in the experimental structure that are also present in
		the predicted structure.*/
		
		/*double totalRealBasePairs = 0;		
		double correctlyPredicted = 0;
		for(int i = 0 ; i < realPairedSites.length ; i++)
		{
			if(realPairedSites[i] >= i + 1 && realPairedSites[i] != 0)
			{
				totalRealBasePairs++;
				if(realPairedSites[i] == predictedPairedSites[i])
				{
					correctlyPredicted++;
				}
			}
		}
		

		//System.out.println(correctlyPredicted);
		//System.out.println(totalRealBasePairs);
		
		//return correctlyPredicted / totalRealBasePairs;
		
		double [] ret = getValues(realPairedSites, predictedPairedSites);
		return ret[0] / (ret[0]+ret[3]);*/
		
		double count = 0;
		double total = 0;
		for(int i = 0 ; i < realPairedSites.length ; i++)
		{
			if(realPairedSites[i] != 0)
			{
				total++;
				if(realPairedSites[i] == predictedPairedSites[i])
				{
					count++;
				}
			}
		}
		
		return count / total;
	}
	
	/**
	 * Given an array of paired sites corresponding to the real structure 
	 * and an array corresponding to the predicted structure returns the PPV.
	 * @param realPairedSites
	 * @param predictedPairedSites
	 */
	public static double calculatePPV (int [] realPairedSites, int [] predictedPairedSites)
	{	
		/* The PPV is the percentage of base
		pairs in the predicted structure that are in the experimental structure.*/
		
		/*
		double correctlyPredicted = 0; // true positives
		double incorrectlyPredicted = 0; // false positives
		for(int i = 0 ; i < realPairedSites.length ; i++)
		{
			if(realPairedSites[i] >= i + 1 && realPairedSites[i] != 0 && realPairedSites[i] == predictedPairedSites[i])
			{
				correctlyPredicted++;
			}
			else
			if(predictedPairedSites[i] >= i + 1 && predictedPairedSites[i] != 0 && realPairedSites[i] != predictedPairedSites[i])
			{
				incorrectlyPredicted++;
			}
		}
		
		//return correctlyPredicted / (correctlyPredicted + incorrectlyPredicted);

		double [] ret = getValues(realPairedSites, predictedPairedSites);
		return ret[0] / (ret[0]+ret[2]);*/
		double count = 0;
		double total = 0;
		for(int i = 0 ; i < realPairedSites.length ; i++)
		{
			if(predictedPairedSites[i] != 0)
			{
				total++;
				if(predictedPairedSites[i] == realPairedSites[i])
				{
					count++;
				}
			}
		}
		
		return count / total;
	}
	
	/**
	 * Given an array of paired sites corresponding to the real structure 
	 * and an array corresponding to the predicted structure returns the F-score.
	 * @param realPairedSites
	 * @param predictedPairedSites
	 */
	public static double calculateFScore (int [] realPairedSites, int [] predictedPairedSites)
	{
		double sensitivity = calculateSensitivity(realPairedSites, predictedPairedSites);
		double ppv = calculatePPV(realPairedSites, predictedPairedSites);
		
		return (2 * sensitivity * ppv)/(sensitivity+ppv);
	}
	
	public static double [] getValues(int [] realPairedSites, int [] predictedPairedSites)
	{
		double TP = 0; // true positives
		double TN = 0;
		double FP = 0; // false positives
		double FN  = 0;
		for(int i = 0 ; i < realPairedSites.length ; i++)
		{
			if(i + 1 < realPairedSites[i] && realPairedSites[i] == predictedPairedSites[i])
			{
				TP++;
			}
			else
			if(i + 1 < realPairedSites[i] && realPairedSites[i] != 0 && predictedPairedSites[i] != realPairedSites[i])
			{
				FN++;
			}
			else
			if(predictedPairedSites[i] != 0 && realPairedSites[i] == 0)
			{
				FP++;
			}
			else
			if(realPairedSites[i] == 0 && predictedPairedSites[i] == 0)
			{
				TN++;
			}
		}
		double [] ret = {TP, TN, FP, FN};
		return ret;		
	}
	
	public static void printValues(int [] realPairedSites, int [] predictedPairedSites)
	{
		double [] ret = getValues(realPairedSites, predictedPairedSites);
		System.out.println("TP="+ret[0]+"  TN="+ret[1]+"   FP="+ret[2]+"   FN="+ret[3]);
	}
	

    
    public static StatAlignResult loadStatAlignResultFile(File resultFile)
    {
    	try
    	{
	    	 BufferedReader buffer = new BufferedReader(new FileReader(resultFile));
	    	 StatAlignResult result = new StatAlignResult();
	    	 result.sequence = buffer.readLine();
	    	 result.dbnStructure = buffer.readLine();
	    	 result.pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketString(result.dbnStructure, '(', ')');
	    	 buffer.close();
	    	 return result;
    	}
    	catch(IOException ex)
    	{
    		ex.printStackTrace();
    	}
    	return null;
    }
    
    static class Scores
    {
    	int [] pairedExperimental;
    	int [] pairedPredicted;
    	double fsc;
    	double sen;
    	double ppv;
    	
    	public static Scores getScores(int [] pairedExperimental, int [] pairedPredicted)
    	{
    		Scores scores = new Scores();
    		scores.pairedExperimental = pairedExperimental;
    		scores.pairedPredicted = pairedPredicted;
    		scores.fsc = Benchmarks.calculateFScore(pairedExperimental, pairedPredicted);
    		scores.sen = Benchmarks.calculateSensitivity(pairedExperimental, pairedPredicted);
    		scores.ppv = Benchmarks.calculatePPV(pairedExperimental, pairedPredicted);
    		
    		if(Double.isNaN(scores.fsc))
    		{
    			scores.fsc = 0;
    		}
    		if(Double.isNaN(scores.ppv))
    		{
    			scores.ppv = 0;
    		}
    		if(Double.isNaN(scores.sen))
    		{
    			scores.sen = 0;
    		}
    		return scores;
    	}
    }
    
    public static double IQR(ArrayList<Double> values)
    {
    	return getValue(values, 0.75) - getValue(values, 0.25);
    }
    
    public static double getMedian(ArrayList<Double> values)
    {
    	ArrayList<Double> sortedValues = (ArrayList<Double>) values.clone();
    	Collections.sort(sortedValues);
    	
    	double length = ((double)values.size()-1)/2;
    	int floor = (int) Math.floor(length);
    	int ceil = (int) Math.ceil(length);
    	
    	return (sortedValues.get(floor)+sortedValues.get(ceil))/2;
    }
    
    public static double getValue(ArrayList<Double> values, double percentile)
    {
    	ArrayList<Double> sortedValues = (ArrayList<Double>) values.clone();
    	Collections.sort(sortedValues);
    	
    	int lower = (int)(percentile * ((double)sortedValues.size()-1));
    	int upper = (int)Math.ceil(percentile * ((double)sortedValues.size()-1));
    	
    	return (sortedValues.get(lower)+sortedValues.get(upper))/2;
    }
    
    public double percentile(ArrayList<Double> values, double x)
    {
    	ArrayList<Double> sortedValues = (ArrayList<Double>) values.clone();
    	Collections.sort(sortedValues);
    	double minIndex = 0;
    	for(int i = 0 ; i < sortedValues.size() ; i++)
    	{
    		if(x >= sortedValues.get(i))
    		{
    			minIndex = i+1;
    			break;
    		}
    	}
    	double maxIndex = sortedValues.size();
    	for(int i = sortedValues.size() - 1 ; i >=  0 ; i--)
    	{
    		if(x <= sortedValues.get(i))
    		{
    			maxIndex = i+1;
    			break;
    		}
    	}
    	
    	double pos = (minIndex+maxIndex)/2;    	
    	return pos / ((double)(sortedValues.size()));
    }
}
