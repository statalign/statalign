package statalign.postprocess.plugins.benchmarks;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import statalign.postprocess.plugins.PPFold;
import statalign.postprocess.utils.Mapping;
import statalign.postprocess.utils.RNAFoldingTools;


public class Benchmarks 
{
	public static void main(String[] args) {
		new Benchmarks().performDistanceBenchmarks();
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
	
	public static void automatedTest2()
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
					//System.out.println(structures);
				}
			}
		}
	}
	
	public void performDistanceBenchmarks()
	{
		File distanceFile = new File("/home/michael/workspace/StatAlign/3seq/dist_scores.txt");
		//String dataDir = "/home/michael/Dropbox/RNA and StatAlign/TestRNAData/";
		
		String dataDir = "/home/michael/Dropbox/RNA and StatAlign/Distance/Datasets2/";
		String resultsDir = "/home/michael/workspace/StatAlign/3seq/";
		try
		{
			BufferedReader buffer = new BufferedReader(new FileReader(distanceFile));
			String textline = null;
			while((textline = buffer.readLine()) != null)
			{
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
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
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
}
