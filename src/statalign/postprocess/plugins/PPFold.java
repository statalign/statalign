package statalign.postprocess.plugins;

import java.awt.BorderLayout;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.base.InputData;
import statalign.base.Mcmc;
import statalign.base.State;
import statalign.base.Utils;
import statalign.distance.Distance;
import statalign.postprocess.Postprocess;
import statalign.postprocess.PostprocessManager;
import statalign.postprocess.gui.PPFoldGUI;
import statalign.postprocess.plugins.benchmarks.Benchmarks;
import statalign.postprocess.plugins.benchmarks.Dataset;
import statalign.postprocess.utils.Mapping;
import statalign.postprocess.utils.RNAFoldingTools;
import statalign.postprocess.utils.RNAalifold;

import com.ppfold.algo.AlignmentData;
import com.ppfold.algo.FuzzyAlignment;
import com.ppfold.algo.NullProgress;
import com.ppfold.algo.Parameters;
import com.ppfold.algo.Progress;
import com.ppfold.algo.ResultBundle;
import com.ppfold.algo.Tree;
import com.ppfold.algo.extradata.ExtraData;
import com.ppfold.main.Alignment;
import com.ppfold.main.AlignmentReader;
import com.ppfold.main.NewickReader;
import com.ppfold.main.PPfoldMain;

public class PPFold extends statalign.postprocess.Postprocess {

	public static void main(String [] args)
	{
		//System.out.println(CalculateMCMCParameters.stepSize());
		
		ArrayList < ArrayList< ArrayList<String> > > allTests = readAllSamples("/home/ingolfur/oxford_workspace/cmdStatAlign/output/alisamples");
		for(int i = 0; i<allTests.size(); ++i){
			ArrayList< ArrayList<String> > currTest= allTests.get(i);
			ArrayList<String> firstSample = currTest.get(0);
			int size = 0;
			int gaps = 0;
			for(String row : firstSample){
				for(int j = 0; j<row.length(); ++j){
					size++;
					if(row.charAt(j) == '-' || row.charAt(j) == '.'){
						gaps++;
					}
				}
			}
			System.out.println((double)gaps/size);
		}
		
	}
	
	
	
	public static ArrayList < ArrayList< ArrayList<String> > > readAllSamples(String path){
		try {
			ArrayList < ArrayList< ArrayList<String> > > allTests = new ArrayList<ArrayList<ArrayList<String>>>();
			//File dir = new File("/home/ingolfur/oxford_workspace/cmdStatAlign/output/alisamples");
			File dir = new File(path);
			File[] sortedDir = dir.listFiles();
			Arrays.sort(sortedDir);
			ArrayList<String> Alignment = new ArrayList<String>();
			ArrayList< ArrayList<String> > Samples= new ArrayList<ArrayList<String>>();


			for(File i : sortedDir){
				Samples = new ArrayList<ArrayList<String>>();
				BufferedReader br = new BufferedReader(new FileReader (i.getAbsoluteFile()));
				String line = br.readLine();
				do{
					if(line.compareTo("%") == 0){
						Samples.add(Alignment);
						Alignment = new ArrayList<String>();
					}else {
						Alignment.add(line);
					}



				}while( (line = br.readLine()) != null);
				allTests.add(Samples);
				br.close();
			}
		return allTests;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	

	// variables from ppfoldmain

	static private Progress progress = NullProgress.INSTANCE;

	public String title;
	
	JPanel pan = new JPanel(new BorderLayout());
	PPFoldGUI gui;
	// private boolean sampling = true;

	CurrentAlignment curAlig;
	MpdAlignment mpdAlignment;

	ColumnNetwork network;
	Column firstVector, lastVector;
	int sizeOfAlignments;

	int[] firstDescriptor;
	static String t[][];
	String[] sequences;
	String[] viterbialignment;
	int d;

	int seqNo = 0;
	static String refSeqName = "";
	String refSeq;
	public String refSeqGapped;

	double[][] summedBasePairProbMatrix;
	double[] summedSingleBaseProb;
	public float[][] probMatrix;	
	double[][] summedBasePairProbRNAalifold;

	int noSamples;

	double [][] weightedBasePairProb;
	double beta = 10;
	double weightedSum = 0;
	double firstLikelihood = 0;
	double posteriorProbabilityAvg = 0;
	ArrayList<String> tempAlignment;
	ArrayList<Double> logLikelihood;
	public double entropyObs;
	public double entropyExp;
	public double entropySample;
	

	RNAFoldingTools rnaTools = new RNAFoldingTools();
	ArrayList<Double> distanceList;

	String outDir = "output/";
	
	boolean samplingAndAveragingPPfold = false;
	boolean samplingAndAveragingRNAalifold = false;
	boolean fuzzyFolding = false;
	boolean experimental = false;

	public PPFold() {
		new File(outDir).mkdirs();
		screenable = true;
		outputable = true;
		postprocessable = true;
		postprocessWrite = true;
		rnaAssociated = true;
		
		if(!show){
			PostprocessManager.rnaMode = true;	
		}
		
		
		Entropy.allowed = true;
	}

	@Override
	public String getTabName() {
		return "PPFold";
	}

	@Override
	public Icon getIcon() {
		return new ImageIcon(ClassLoader.getSystemResource("icons/MPD.gif"));
	}

	@Override
	public JPanel getJPanel() {
		return pan;
	}

	@Override
	public String getTip() {
		return "Base-pairing matrix of the current consensus structure given by PPFold";
	}
	
	@Override
    public double getTabOrder() {
        return 8.0d;
    }

	@Override
	public void setSampling(boolean enabled) {
		sampling = enabled;
	}

	@Override
	public String[] getDependences() {
		return new String[] { "statalign.postprocess.plugins.CurrentAlignment", "statalign.postprocess.plugins.MpdAlignment"};
	}

	@Override
	public void refToDependences(Postprocess[] plugins) {
		curAlig = (CurrentAlignment) plugins[0];
		mpdAlignment = (MpdAlignment) plugins[1];
	}

	static Comparator<String[]> compStringArr = new Comparator<String[]>() {
		@Override
		public int compare(String[] a1, String[] a2) {
			return a1[0].compareTo(a2[0]);
		}
	};



	List<ExtraData> extradata = new ArrayList<ExtraData>();
	Parameters param = null;

	
	public void selectReferenceSequence(InputData input)
	{
		int maxLength = 0;
		refSeq = input.seqs.sequences.get(0).replaceAll("-", "");
		refSeqName = input.seqs.seqNames.get(0);
		refSeqGapped = input.seqs.sequences.get(0);
		maxLength = refSeq.length();
		for (int i = 0; i < input.seqs.sequences.size(); i++) {
			String seq = input.seqs.sequences.get(i).replaceAll("-", "");
			if (seq.length() > maxLength) {
				maxLength = seq.length();
				refSeq = seq;
				refSeqName = input.seqs.seqNames.get(i);
				refSeqGapped = input.seqs.sequences.get(i);
				seqNo = i;
			}
		}
		d = maxLength;
	}

	@Override
	public void beforeFirstSample(InputData input) {		
		title = input.title + "_seed"+input.pars.seed;
		selectReferenceSequence(input);
		
		System.out.println("IS NULL = " +pluginParameters);
		if(pluginParameters != null)
		{
			System.out.println("Parameters");
			pluginParameters.print();
			String ppfoldParameter = pluginParameters.getParameter("ppfold");
			String rnaAlifoldParameter = pluginParameters.getParameter("rnaalifold");
			String fuzzyParameter = pluginParameters.getParameter("fuzzy");
			
			samplingAndAveragingPPfold = ppfoldParameter != null;
			samplingAndAveragingRNAalifold = rnaAlifoldParameter != null;
			fuzzyFolding = true;
			
			if(samplingAndAveragingRNAalifold)
			{
				System.out.println("Using RNAalifold with parameters: \"" + rnaAlifoldParameter+"\"");
			}
			
			if(pluginParameters.getParameter("experimental") != null)
			{
				samplingAndAveragingPPfold = true;
				samplingAndAveragingRNAalifold = true;
				fuzzyFolding = true;
				experimental = true;
			}

		}
		
		if(show) {
			pan.removeAll();
			JScrollPane scroll = new JScrollPane();
			scroll.getViewport().add(gui = new PPFoldGUI(title, scroll, this));

			if(gui.getPreferredSize().getHeight() > gui.getHeight()) {
				scroll.createVerticalScrollBar();
			}

			if(gui.getPreferredSize().getWidth() > gui.getWidth()) {
				scroll.createHorizontalScrollBar();
			}

			pan.add(scroll, BorderLayout.CENTER);
			pan.getParent().validate();
			
			gui.changeDimension(d);

		}
		
		sizeOfAlignments = (mcmc.tree.vertex.length + 1) / 2;
		noSamples = 0;

		t = new String[sizeOfAlignments][];
		sequences = null;
		viterbialignment = new String[sizeOfAlignments];

		network = new ColumnNetwork();

		firstDescriptor = new int[sizeOfAlignments];
		Arrays.fill(firstDescriptor, -1);
		firstVector = network.add(firstDescriptor);
		lastVector = null;

		viterbialignment = new String[sizeOfAlignments];
		logLikelihood = new ArrayList<Double>();

		// read in PPfold parameter file
		try {
			BufferedReader paramFileReader = null;
			if(param == null)
			{
				File file = new File("res/matrices.in");
				paramFileReader = new BufferedReader(new InputStreamReader(	new FileInputStream(file)));
				param = Parameters.readParam(paramFileReader);				
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		// initialise new dataset
		dataset = new Dataset();
			

		String [] [] inputAlignment = new String[input.seqs.sequences.size()][2];

		for(int k = 0 ; k < inputAlignment.length ; k++)
		{
			inputAlignment[k][0] = input.seqs.seqNames.get(k);
			inputAlignment[k][1] = input.seqs.sequences.get(k);
		}
		Arrays.sort(inputAlignment, compStringArr);

		dataset.inputAlignment = new AlignmentData();
		for(int i = 0 ; i < inputAlignment.length ; i++)
		{
			dataset.inputAlignment.names.add(inputAlignment[i][0]);
			dataset.inputAlignment.sequences.add(inputAlignment[i][1]);
		}		


		if(experimental)
		{
						// perform calculation on reference sample
			try {
				
				ResultBundle refResult = PPfoldMain.fold(progress, input.seqs.sequences, input.seqs.seqNames, null, param, extradata);
				int [] refPairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(refResult.getFinalMatrix());
				dataset.resultBundlePPfold = refResult.getSmallBundle();
				dataset.pairedSitesPPfoldProjected = Benchmarks.projectPairedSites(RNAFoldingTools.getSequenceByName(refSeqName, input.seqs.sequences, input.seqs.seqNames), refPairedSites);
				
				dataset.rnaAlifoldRef = RNAalifold.fold(input.seqs.sequences, input.seqs.seqNames, Postprocess.pluginParameters.getParameter("rnaalifold")).getSmallResult();
				dataset.pairedSitesRNAalifoldRefProjected = Benchmarks.projectPairedSites(RNAFoldingTools.getSequenceByName(refSeqName, input.seqs.sequences, input.seqs.seqNames), dataset.rnaAlifoldRef.pairedSites);
				//dataset.pairedSitesRNAalifoldRefProjected = Benchmarks.projectPairedSites(refSeq, dataset.rnaAlifoldRef.pairedSites);
				
				//RNAFoldingTools.writeToFile(new File(outDir+title+"_ref_structure.txt"), ""+refResult.getPPfoldReliability()+"\n"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites)+"\n", false);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		alignments = new ArrayList<AlignmentData>();
	}
	
	public double calculateMPD(ArrayList<String> mpdSequences, ArrayList<String> mpdNames)
	{
		try
		{
			List<String> sequences = mpdSequences;
			List<String> seqNames = mpdNames;
			String refSeq = sequences.get(0).replaceAll("-", "");
			String refSeqGapped = sequences.get(0);
			int maxLength = refSeq.length();
			for (int i = 0; i < sequences.size(); i++) {
				String seq = sequences.get(i).replaceAll("-", "");
				if (seq.length() > maxLength) {
					maxLength = seq.length();
					refSeq = seq;
					refSeqName = seqNames.get(i);
					refSeqGapped = sequences.get(i);
				}
			}
			

			List<ExtraData> extradata = new ArrayList<ExtraData>();
			
			ResultBundle mpdResult = PPfoldMain.fold(progress, mpdSequences, mpdNames, null, param, extradata);
			float [][] mpdBasePairProb = mpdResult.finalmatrix;			
			
			int [] mpdPairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(mpdBasePairProb);
			//double ppfoldReliability2 = RNAFoldingTools.calculatePPfoldReliabilityScore(mpdPairedSites, RNAFoldingTools.getDoubleMatrix(mpdBasePairProb));
			
			
			double ppfoldReliability = mpdResult.getPPfoldReliability();
			int [] projectedPairedSites = Benchmarks.projectPairedSites(refSeqGapped, mpdPairedSites);
			dataset.resultBundleMPD = mpdResult.getSmallBundle();
			dataset.pairedSitesMPD = projectedPairedSites;
			dataset.ppfoldReliabilityMPD = ppfoldReliability;
			dataset.pairsOnlyReliabilityMPD = mpdResult.getPairsOnlyReliability();
			dataset.pairsOnlyMPDPosteriorWeighted = RNAFoldingTools.calculatePairsOnlyReliabilityScore(mpdPairedSites, RNAFoldingTools.getDoubleMatrix(mpdResult.finalmatrix), dataset.posteriors);
			
			RNAalifoldResult rnaAlifoldMPDResult = RNAalifold.fold(mpdSequences, mpdNames, Postprocess.pluginParameters.getParameter("rnaalifold"));
			dataset.rnaAlifoldMPD = rnaAlifoldMPDResult.getSmallResult();
			dataset.pairedSitesRNAalifoldMPDProjected = Benchmarks.projectPairedSites(refSeqGapped, rnaAlifoldMPDResult.pairedSites);
	
			return ppfoldReliability;
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
		}

		return -1;
	}


	Dataset dataset = null;
	ArrayList<AlignmentData> alignments = new ArrayList<AlignmentData>();
	public AlignmentData al;
	//String p2 = "";
	//double finalEntropyExpReliabilityScore = -1;
	//double finalEntropyObsReliabilityScore = -1;
	
	@Override	
	public void newSample(State state, int no, int total) {
		for (int i = 0; i < t.length; i++) {
			t[i] = curAlig.leafAlignment[i].split("\t");
		}
		Arrays.sort(t, compStringArr);	
		
		if(experimental)
		{
			if(no == 0)
			{
				dataset.title = title;
				dataset.randomSeed = mcmc.mcmcpars.seed;
				dataset.burnIn = mcmc.mcmcpars.burnIn;
				dataset.mcmcSteps = mcmc.mcmcpars.cycles;
				dataset.samplingRate = mcmc.mcmcpars.sampRate;
			}
		}
		
		if (sequences == null) {
			sequences = new String[sizeOfAlignments];		
			int len = t[0][1].length();
			for (int i = 0; i < sizeOfAlignments; i++) {
				sequences[i] = "";
				for (int j = 0; j < len; j++) {
					if (t[i][1].charAt(j) != '-') {
						sequences[i] += t[i][1].charAt(j);
					}
				}
			}
		}


		
		if(no == 0)
		{
			if(samplingAndAveragingPPfold)
			{
				summedBasePairProbMatrix = new double[d][d];
				weightedBasePairProb = new double[d][d];
				summedSingleBaseProb = new double[d];
			}
			if(samplingAndAveragingRNAalifold)
			{
				summedBasePairProbRNAalifold = new double[d][d];			
			}
		}

		
		boolean cont = true;
		if (cont) {
			/*// initialise matrices
			if (no == 0) {
			
				weightedBasePairProb = new double[d][d];
				for (int i = 0; i < d; ++i) {
					for (int j = 0; j < d; ++j) {
						summedBasePairProbMatrix[i][j] = 0;
						weightedBasePairProb[i][j] = 0;
						summedBasePairProbRNAalifold[i][j] = 0;
					}
				}
				summedSingleBaseProb = new double[d];
			}*/

			try {
				List<String> lines = new ArrayList<String>();
				for (int i = 0; i < sequences.length; ++i) {
					lines.add(">" + t[i][0].trim());
					lines.add(t[i][1]);
				}

				Alignment align = AlignmentReader.readAlignmentFromStringList(lines);

			

				Tree tree = getPPfoldTree(mcmc);
				
				if(samplingAndAveragingPPfold)
				{
					ResultBundle sampleResult = PPfoldMain.fold(progress, align.getSequences(),	align.getNames(), tree, param, extradata);
					entropySample = sampleResult.entropyVal;
					float[][] basePairProb = sampleResult.finalmatrix;
					float[] singleBaseProb = new float[basePairProb.length];
					for (int x = 0; x < basePairProb.length; x++) {
						singleBaseProb[x] = 1;
						for (int y = 0; y < basePairProb[0].length; y++) {
							singleBaseProb[x] -= basePairProb[x][y];
						}
					}


					ArrayList<String> sequences = new ArrayList<String>();
					for (int k = 0; k < t.length; k++) {
						sequences.add(t[k][1]);
					}
					float[][] projectSample = Mapping.projectMatrix(PPFold.getSequenceByName(t, refSeqName), basePairProb, '-');
					float[] projectSingleBaseProb = Mapping.projectarray(PPFold.getSequenceByName(t, refSeqName),singleBaseProb, '-');

					// normalise projected matrix
					for(int x = 0 ; x < projectSample.length ; x++)
					{
						double rowMatrSum = 0;
						for(int y = 0 ; y < projectSample[0].length ; y++)
						{
							rowMatrSum += projectSample[x][y];
						}
						double factor = rowMatrSum + projectSingleBaseProb[x];
						factor = 1;
						for(int y = 0 ; y < projectSample[0].length ; y++)
						{
							projectSample[x][y] = (float)(projectSample[x][y] / factor);
						}
						projectSingleBaseProb[x] /= factor;
					}

					
					double alignmentLogLikelihood = state.logLike;
					
					
					probMatrix = new float[d][d];
					//float [][] rnaAlifoldProbMatrix = new float[d][d];
					double weight = Math.pow(firstLikelihood / alignmentLogLikelihood, beta);

					float [] singleMatrix = new float[d];
					for (int i = 0; i < d; ++i) {
						summedSingleBaseProb[i] += projectSingleBaseProb[i];
						singleMatrix[i] = (float)summedSingleBaseProb[i] / (noSamples+1);
						for (int j = 0; j < d; ++j) {
							summedBasePairProbMatrix[i][j] += projectSample[i][j];
							probMatrix[i][j] = (float)(summedBasePairProbMatrix[i][j]/(noSamples+1));
							weightedBasePairProb[i][j] += projectSample[i][j]*weight;
						}
					}
					weightedSum += weight;

					
					Structure.updateBasePairMatrix(probMatrix);
					Structure.updateSingleMatrix(singleMatrix);
					PPfoldMain.setfoldingfinished(true);
					
					if(gui != null)
					{
						System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(probMatrix)));
						System.out.println("PAINTING");
						gui.changeDimension(d*PPFoldGUI.OFFSET);
						gui.setMatrix(probMatrix);
						gui.repaint();
					}
					
					if(experimental)
					{
						dataset.sampledStructures.add(sampleResult.getSmallBundle());
						int [] samplePairedSitesProjected =  rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(projectSample);						
						dataset.pairedSitesProjectedSamples.add(samplePairedSitesProjected);
					}
				}
				
				if(samplingAndAveragingRNAalifold)
				{
					RNAalifoldResult rnaAlifoldSampleResult = RNAalifold.fold(align.getSequences(), align.getNames(), Postprocess.pluginParameters.getParameter("rnaalifold"));
					double [][] rnaAlifoldMatrixSample = rnaAlifoldSampleResult.matrix;

					ArrayList<String> sequences = new ArrayList<String>();
					for (int k = 0; k < t.length; k++) {
						sequences.add(t[k][1]);
					}


					float[][] rnaAlifoldProjectedSample = Mapping.projectMatrix(PPFold.getSequenceByName(t, refSeqName), RNAFoldingTools.getFloatMatrix(rnaAlifoldMatrixSample), '-');
					for (int i = 0; i < d; ++i) {
						for (int j = 0; j < d; ++j) {
							summedBasePairProbRNAalifold[i][j] += rnaAlifoldProjectedSample[i][j];
						}
					}
					
					/*
					Structure.updateMatrix(probMatrix);					
					PPfoldMain.setfoldingfinished(true);
					
					if(show)
					{
						gui.changeDimension(d*PPFoldGUI.OFFSET);
						gui.setMatrix(probMatrix);
						gui.repaint();
					}*/
					
					if(experimental)
					{
						dataset.pairedSitesProjectedRnaAlifoldSamples.add(rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(rnaAlifoldProjectedSample));
					}					
				}				

			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			noSamples += 1;
		}

		// save alignment sample information
		if(no == 0)
		{
			

			AlignmentData al =  new AlignmentData();
			for(int k = 0 ; k < t.length ; k++)
			{
				al.sequences.add(t[k][1]);
				al.names.add(t[k][0]);
			}
			alignments.add(al);

			/*
			appendAlignment("reference", inputAlignment, new File(outDir+mpdAlignment.input.title+".samples"), false);
			appendAlignment(no+"", t, new File(outDir+mpdAlignment.input.title+".samples"), true);
			RNAFoldingTools.writeToFile(new File(outDir +title+"_Distance.txt"),"BurnIn "+new Integer(mcmc.mcmcpars.burnIn).toString() + "\t" +
					"cycles "+new Integer(mcmc.mcmcpars.cycles).toString() + "\t" +
					"sampRate "+new Integer(mcmc.mcmcpars.sampRate).toString() + "\t" + 
					"Filename "+ mpdAlignment.input.title + "\t", false);*/


			tempAlignment = new ArrayList<String>();
			for(int k = 0 ; k < t.length ; k++)
			{
				tempAlignment.add(t[k][1]);
			}
			
			//& appendAlignment("reference", inputAlignment, new File(outDir+title+".samples"), false);
			//& appendAlignment(no+"", t, new File(outDir+title+".samples"), true);
			
			//& RNAFoldingTools.writeToFile(new File(outDir+title+"_entropy_fuzzy_obs.txt"),"", false);
			//& RNAFoldingTools.writeToFile(new File(outDir+title+"_entropy_fuzzy_exp.txt"),"", false);
			//& RNAFoldingTools.writeToFile(new File(outDir+title+"_entropy_samples.txt"),"", false);
			//& RNAFoldingTools.writeToFile(new File(outDir+title+"_entropy_mpd.txt"),"", false);

			





		}

		if(experimental)
		{
			dataset.logLikelihoods.add(mcmc.mcmcStep.newLogLike);
		}
		
		if(fuzzyFolding)
		{
			al =  new AlignmentData();
			for(int k = 0 ; k < t.length ; k++)
			{
				al.sequences.add(t[k][1]);
				al.names.add(t[k][0]);
			}
			alignments.add(al);

			try
			{
				FuzzyAlignment fuzzyAlignment2 = FuzzyAlignment.getFuzzyAlignmentAndProject(alignments, refSeqName);
				ResultBundle fuzzyResultExp2 = PPfoldMain.foldFuzzyAlignment(progress, fuzzyAlignment2, null, param, extradata, true);
				ResultBundle fuzzyResultObs2 = PPfoldMain.foldFuzzyAlignment(progress, fuzzyAlignment2, null, param, extradata, false);
				entropyObs = fuzzyResultObs2.entropyVal;
				entropyExp = fuzzyResultExp2.entropyVal;
			}
			catch(Exception ex)
			{
				ex.printStackTrace();
			}
			
			
			
			if(experimental)
			{
				double dist = Distance.AMA(tempAlignment, al.sequences);
				RNAFoldingTools.writeToFile(new File(outDir +title+"_Distance.txt"),new Double(dist).toString(), true);
				
				FuzzyAlignment fuzzyAlignment = FuzzyAlignment.getFuzzyAlignmentAndProject(alignments, refSeqName);
				
				List<ExtraData> extradata = new ArrayList<ExtraData>();
				try {
					Tree tree = null;
					//Tree tree = getPPfoldTree(mcmc);				
					//getPPfoldTree(mcmc).print();
					
					ResultBundle fuzzyResultExp = PPfoldMain.foldFuzzyAlignment(progress, fuzzyAlignment, tree, param, extradata, true);
					ResultBundle fuzzyResultObs = PPfoldMain.foldFuzzyAlignment(progress, fuzzyAlignment, tree, param, extradata, false);
					double finalEntropyExpReliabilityScore = fuzzyResultExp.getPPfoldReliability();
					double finalEntropyObsReliabilityScore = fuzzyResultObs.getPPfoldReliability();
					int [] fuzzyPairedSitesExp = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(fuzzyResultExp.finalmatrix);
					int [] fuzzyPairedSitesObs = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(fuzzyResultObs.finalmatrix);
					dataset.resultBundleEntropyExp = fuzzyResultExp.getSmallBundle();
					dataset.resultBundleEntropyObs = fuzzyResultObs.getSmallBundle();
					dataset.pairedSitesEntropyExp = fuzzyPairedSitesExp;
					dataset.pairedSitesEntropyObs = fuzzyPairedSitesObs;				
					dataset.cumulativeFuzzyAlignment.add(fuzzyAlignment);
					dataset.sampledAlignments.add(al);
					dataset.cumulativeFuzzyExpResults.add(fuzzyResultExp.getSmallBundle());
					dataset.cumulativeFuzzyObsResults.add(fuzzyResultObs.getSmallBundle());
					
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}

	}
	int stepCounter = 0;
	
	public void computeFuzzyAlignment()
	{
		FuzzyAlignment fuzzyAlignment = FuzzyAlignment.getFuzzyAlignmentAndProject(alignments, refSeqName);
		System.out.println("Fuzzy alignment size"+fuzzyAlignment.columns.size());
		List<ExtraData> extradata = new ArrayList<ExtraData>();
		try {
			Tree tree = null;
			ResultBundle fuzzyResultObs = PPfoldMain.foldFuzzyAlignment(progress, fuzzyAlignment, tree, param, extradata, false);
			System.out.println("Fuzzy matrix size"+fuzzyResultObs.finalmatrix[0].length);
			int [] fuzzyPairedSitesObs = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(fuzzyResultObs.finalmatrix);			
			char [] structure = fuzzyResultObs.getStructure();
			
			RNAFoldingTools.saveCtFile(new File(title+".fuzzy.ct"), fuzzyPairedSitesObs, title, refSeq);
			RNAFoldingTools.saveDotBracketFile(new File(title+".fuzzy.dbn"), fuzzyPairedSitesObs, title, refSeq);
			
			//if(experimental)
			//{
				//dataset.pp
			dataset.pairedSitesEntropyObs = fuzzyPairedSitesObs;
			dataset.ppfoldReliabilityEntropyObs = RNAFoldingTools.calculatePPfoldReliabilityScore(fuzzyPairedSitesObs, RNAFoldingTools.getDoubleMatrix(fuzzyResultObs.finalmatrix)); 
			dataset.pairsOnlyReliabilityEntropyObs = RNAFoldingTools.calculatePairsOnlyReliabilityScore(fuzzyPairedSitesObs, RNAFoldingTools.getDoubleMatrix(fuzzyResultObs.finalmatrix));
			dataset.pairsOnlyReliabilityEntropyObsPosteriorWeighted = RNAFoldingTools.calculatePairsOnlyReliabilityScore(fuzzyPairedSitesObs, RNAFoldingTools.getDoubleMatrix(fuzzyResultObs.finalmatrix), dataset.posteriors);
			dataset.fuzzyAlignmentObsEntropyVal = fuzzyResultObs.entropyVal;
			dataset.fuzzyAlignmentObsEntropyPerc = fuzzyResultObs.entropyPercOfMax;
			dataset.fuzzyAlignmentObsEntropyMax = fuzzyResultObs.entropyMax;
				
			//}
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}		
	
	}


	@Override
	public void afterLastSample() {
		// calculate posterior avg
		posteriorProbabilityAvg = 0;
		for(int i = 0 ; i < mpdAlignment.decoding.length ; i++)
		{
			posteriorProbabilityAvg += mpdAlignment.decoding[i];
			dataset.posteriors.add(mpdAlignment.decoding[i]);
		}

		posteriorProbabilityAvg /= (mpdAlignment.decoding.length);
		if(experimental)
		{
			dataset.posteriorsAverage = posteriorProbabilityAvg;
		}

		
		
		if(samplingAndAveragingRNAalifold)
		{
			double[][] doubleSummedArrayRNAalifold = new double[d][d];
			for (int i = 0; i < d; ++i) {
				for (int j = 0; j < d; ++j) {
					doubleSummedArrayRNAalifold[i][j] = summedBasePairProbRNAalifold[i][j] / noSamples;
				}
			}

			int [] pairedSitesRNAalifold = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(doubleSummedArrayRNAalifold);
			dataset.pairedSitesRNAalifold = pairedSitesRNAalifold;
			dataset.ppfoldReliabilityScoreRNAalifold = RNAFoldingTools.calculatePPfoldReliabilityScore(dataset.pairedSitesRNAalifold, doubleSummedArrayRNAalifold);
			dataset.pairsOnlyReliabilityScoreRNAalifold =  RNAFoldingTools.calculatePairsOnlyReliabilityScore(dataset.pairedSitesRNAalifold, doubleSummedArrayRNAalifold);		
			
			RNAFoldingTools.saveCtFile(new File(title+".rnaalifold.ct"), pairedSitesRNAalifold, title, refSeq);
			RNAFoldingTools.saveDotBracketFile(new File(title+".rnaalifold.dbn"), pairedSitesRNAalifold, title, refSeq);
		}
		
		
		if(samplingAndAveragingPPfold)
		{
			// average matrices, important for posterior decoding
			double[][] doubleSummedArray = new double[d][d];			
			double[] doubleSingleBaseProb = new double[d];
			for (int i = 0; i < d; ++i) {
				doubleSingleBaseProb[i] = summedSingleBaseProb[i] / noSamples;
				for (int j = 0; j < d; ++j) {
					doubleSummedArray[i][j] = summedBasePairProbMatrix[i][j] / noSamples;			
				}
			}
			int[] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(doubleSummedArray);
			double statalignPpfoldReliablityScore = RNAFoldingTools.calculatePPfoldReliabilityScore(pairedSites, doubleSummedArray);
			RNAFoldingTools.saveCtFile(new File(title+".ppfold.ct"), pairedSites, title, refSeq);
			RNAFoldingTools.saveDotBracketFile(new File(title+".ppfold.dbn"), pairedSites, title, refSeq);
			
			//if(experimental)
			//{
				for (int i = 0; i < d; ++i) {
					for (int j = 0; j < d; ++j) {
						weightedBasePairProb[i][j] = weightedBasePairProb[i][j] / noSamples;			
					}
				}
				
				dataset.pairedSitesWeighted = RNAFoldingTools.getPosteriorDecodingConsensusStructure(weightedBasePairProb);
				dataset.ppfoldReliabilityScoreSamplingAndAveragingWeighted = RNAFoldingTools.calculatePPfoldReliabilityScore(dataset.pairedSitesWeighted, weightedBasePairProb);
				dataset.pairsOnlyReliabilityScoreSamplingAndAveragingWeighted = RNAFoldingTools.calculatePairsOnlyReliabilityScore(dataset.pairedSitesWeighted, weightedBasePairProb);				
				
				dataset.pairedSites = pairedSites;
				dataset.ppfoldReliabilityScoreSamplingAndAveraging = statalignPpfoldReliablityScore;
				dataset.pairsOnlyReliabilityScoreSamplingAndAveraging = RNAFoldingTools.calculatePairsOnlyReliabilityScore(pairedSites, doubleSummedArray);
				dataset.pairsOnlyReliabilityScoreSamplingAndAveragingPosteriorWeighted = RNAFoldingTools.calculatePairsOnlyReliabilityScore(pairedSites, doubleSummedArray, dataset.posteriors);
				dataset.pairedSitesRefSeq = refSeqGapped;
				
				
				for (int i = 0; i < d; ++i) {
					for (int j = 0; j < d; ++j) {
						weightedBasePairProb[i][j] /= weightedSum;
					}
				}
		
				//RNAFoldingTools.writeMatrix(weightedBasePairProb, new File(outDir+"bp_log.matrix"));
				int [] pairedSitesWeighted= rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(weightedBasePairProb);
				double statalignWeightedPpfoldReliablityScore = RNAFoldingTools.calculatePPfoldReliabilityScore(pairedSites, weightedBasePairProb);
				dataset.pairedSitesWeighted = pairedSitesWeighted;
				dataset.ppfoldReliabilityScoreSamplingAndAveragingWeighted = statalignWeightedPpfoldReliablityScore;
				dataset.pairsOnlyReliabilityScoreSamplingAndAveragingWeighted = RNAFoldingTools.calculatePairsOnlyReliabilityScore(pairedSitesWeighted, weightedBasePairProb);
				//& saveResult(refSeqGapped, pairedSitesWeighted, weightedBasePairProb, RNAFoldingTools.getSingleBaseProb(weightedBasePairProb), new File(outDir+title+".dat.res.weighted"));
						
				for(int i = 0 ; i < mpdAlignment.alignment.length && i < dataset.inputAlignment.names.size() ; i++)
				{					
					dataset.mpdAlignment.names.add(dataset.inputAlignment.names.get(i));
					dataset.mpdAlignment.sequences.add(mpdAlignment.alignment[i]);
				}
				dataset.ppfoldReliabilityScoreSamplingAndAveraging = statalignPpfoldReliablityScore;
				
				double ppfoldReliablityScore = RNAFoldingTools.calculatePPfoldReliabilityScore(pairedSites, doubleSummedArray);
				double proxySim = getProxySimilarityFromAvgPosterior(posteriorProbabilityAvg);
				ArrayList<String> mpdSequences = new ArrayList<String>();
				ArrayList<String> mpdNames = new ArrayList<String>();
				//String [] Utils.alignmentTransformation(mpdAlignment.alignment, "Fasta", mpdAlignment.input);
				for(int i = 0 ; i < mpdAlignment.alignment.length ; i++)
				{
					mpdNames.add(mpdAlignment.alignment[i].split("(\\s)+")[0]);
					mpdSequences.add(mpdAlignment.alignment[i].split("(\\s)+")[1]);
				}			
				dataset.mpdVsInputSim = Distance.AMA(mpdSequences, dataset.inputAlignment.sequences);
		
				
				calculateMPD(mpdSequences, mpdNames);
				
				//System.out.println("Proxy distance = " + Distance.amaScoreToMultiDistance(mpdSequences,proxySim));
				double improvedReliabilityScore1 = posteriorProbabilityAvg*statalignPpfoldReliablityScore;		
				//System.out.println("Improved reliability score 1 = " + improvedReliabilityScore1);
				double improvedReliabilityScore2 = proxySim*statalignPpfoldReliablityScore;
			//}
		}	
		
		if(samplingAndAveragingPPfold && samplingAndAveragingRNAalifold)
		{
			double [][] combinedBasePairProb = new double[d][d];
			for(int i = 0 ; i < d ; i++)
			{
				for(int j = 0 ; j < d ; j++)
				{
					combinedBasePairProb[i][j] = (summedBasePairProbMatrix[i][j]+summedBasePairProbRNAalifold[i][j])/2;
				}
			}
			
			int [] combinedPairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(combinedBasePairProb);
			
			dataset.pairedSitesCombined = combinedPairedSites;
			dataset.ppfoldReliabilityScoreCombined = RNAFoldingTools.calculatePPfoldReliabilityScore(combinedPairedSites, combinedBasePairProb);
			dataset.pairsOnlyReliabilityScoreCombined = RNAFoldingTools.calculatePairsOnlyReliabilityScore(combinedPairedSites, combinedBasePairProb);			
		}

		if(fuzzyFolding)
		{		
			computeFuzzyAlignment();
		}
		
		if(experimental)
		{
			dataset.saveDatasetResult(new File(outDir+title+".serialized"));
		}
		
		DecimalFormat df = new DecimalFormat("0.000");
		
		try
		{
			File outFile = new File(title+".info");
			System.out.println("XXX"+outFile);
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile));
			buffer.write(">method=sampling_and_averaging_ppfold\ttitle="+dataset.title+"\tlength="+dataset.pairedSites.length+"\treliability_score="+df.format(dataset.pairsOnlyReliabilityScoreSamplingAndAveraging)+"\n");
			buffer.write(refSeqGapped.replaceAll("-", "")+"\n");
			buffer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(dataset.pairedSites)+"\n");
			buffer.write(">method=sampling_and_averaging_rnaalifold\ttitle="+dataset.title+"\tlength="+dataset.pairedSitesRNAalifold.length+"\treliability_score="+df.format(dataset.pairsOnlyReliabilityScoreRNAalifold)+"\n");
			buffer.write(refSeqGapped.replaceAll("-", "")+"\n");
			buffer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(dataset.pairedSitesRNAalifold)+"\n");
			buffer.write(">method=fuzzy_alignment_ppfold\ttitle="+dataset.title+"\tlength="+dataset.pairedSitesEntropyObs.length+"\treliability_score="+df.format(dataset.pairsOnlyReliabilityEntropyObs)+"\tentropy="+dataset.fuzzyAlignmentObsEntropyVal+"\tperc_of_max_entropy="+dataset.fuzzyAlignmentObsEntropyPerc+"\t"+"\tmax_entropy="+dataset.fuzzyAlignmentObsEntropyMax+"\t"+"\n");
			buffer.write(refSeqGapped.replaceAll("-", "")+"\n");
			buffer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(dataset.pairedSitesEntropyObs)+"\n");
			buffer.close();
			
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}

	public void saveScores(File outFile, double ppfoldReliabilityStatAlign, double ppfoldReliabilityMPD)
	{
		boolean writeHeaders = !outFile.exists();
		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile, !writeHeaders));
			if(writeHeaders)
			{
				buffer.write("Dataset\tPosterior avg.\tProxy similarity\tProxy distance\tPPfold reliability (StatAlign)\tPPfold reliability(MPD)\tBurn-in\tCycles\tSample rate\n");
			}

			double proxySimilarity = getProxySimilarityFromAvgPosterior(posteriorProbabilityAvg);
			ArrayList<String> sequences = new ArrayList<String>();
			for(int i = 0 ; i < mpdAlignment.alignment.length ; i++)
			{
				sequences.add(mpdAlignment.alignment[i]);
			}
			double proxyDistance = Distance.amaScoreToDistance(sequences,proxySimilarity);
			buffer.write(title+"\t");
			buffer.write(posteriorProbabilityAvg+"\t");
			buffer.write(proxySimilarity+"\t");
			buffer.write(proxyDistance+"\t");
			buffer.write(ppfoldReliabilityStatAlign+"\t");
			buffer.write(ppfoldReliabilityMPD+"\t");
			buffer.write(mcmc.mcmcpars.burnIn+"\t");
			buffer.write(mcmc.mcmcpars.cycles+"\t");
			buffer.write(mcmc.mcmcpars.sampRate+"\t");
			buffer.write("\n");
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}

	public static double getProxySimilarityFromAvgPosterior(double avgPosterior)
	{
		return Math.min(1, avgPosterior*0.8084689957 + 0.2011076631);
	}

	public static void appendAlignment(String label, String [][] alignment, File outFile, boolean append)
	{
		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile, append));
			buffer.write("%"+label+"\n");
			for (int i = 0; i < alignment.length; i++) {
				buffer.write(">"+alignment[i][0] + "\n");
				buffer.write(alignment[i][1] + "\n");
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}

	public static void appendAlignment(String label, String [] alignment, File outFile, boolean append, InputData input)
	{
		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile, append));
			buffer.write("%"+label+"\n");
			String [] fastaAlignment = Utils.alignmentTransformation(alignment,"Fasta", input);
			for (int i = 0; i < fastaAlignment.length; i++) {
				buffer.write(fastaAlignment[i] + "\n");
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}

	public static String getSequenceByName(String[][] sequences, String name) {
		for (int i = 0; i < sequences.length; i++) {
			if (sequences[i] != null && sequences[i][0].equals(name)) {
				return sequences[i][1];
			}
		}
		return null;
	}

	public static void saveResult(String sequence, int[] pairedSites,
			double[][] basePairProb, double[] singleBaseProb, File outFile) {
		DecimalFormat df = new DecimalFormat("0.0000");
		try {
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile));
			buffer.write(sequence + "\n");
			buffer.write(RNAFoldingTools
					.getDotBracketStringFromPairedSites(pairedSites) + "\n");
			for (int k = 0; k < pairedSites.length; k++) {
				if (pairedSites[k] == 0) {
					buffer.write(RNAFoldingTools.pad((k + 1) + "", 4)
							+ "\t"
							+ RNAFoldingTools.pad(pairedSites[k] + "", 7)
							+ RNAFoldingTools.pad("-", 6)
							+ "\t"
							+ RNAFoldingTools.pad(df.format(singleBaseProb[k])
									+ "", 6) + "\n");
				} else {
					buffer.write(RNAFoldingTools.pad((k + 1) + "", 4)
							+ "\t"
							+ RNAFoldingTools.pad(pairedSites[k] + "", 7)
							+ RNAFoldingTools.pad(
									df.format(basePairProb[k][pairedSites[k] - 1]),
									6)
									+ "\t"
									+ RNAFoldingTools.pad(df.format(singleBaseProb[k])
											+ "", 6) + "\n");
				}
			}

			buffer.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		}
	}

	public static com.ppfold.algo.Tree getPPfoldTree(Mcmc mcmc) {
		try {
			return NewickReader.parse(mcmc.tree.printedTree());
		} catch (Exception ex) {
			ex.printStackTrace();
		}

		return null;
	}

	public static String[][] getSequences() {
		return t;
	}

	public static String getRefName() {
		return refSeqName;
	}

	public static void appendFolds(File file, String name, String alignedSequence, int [] pairedSites, int [] projectedPairedSites, boolean append)
	{

		try {
			BufferedWriter buffer = new BufferedWriter(new FileWriter(file, append));
			buffer.write(">"+name+"\n");
			buffer.write(alignedSequence+"\n");
			buffer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites)+"\n");
			buffer.write(alignedSequence.replaceAll("-", "")+"\n");
			buffer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(projectedPairedSites)+"\n");
			buffer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private boolean getTheNextSample(ArrayList<Double> distanceList, int sampleNumber){
		if(sampleNumber > 40 && isThereMajorIncrease(distanceList)){
			return true;
		}
		else{
			return false;
		}
	}

	private boolean isThereMajorIncrease(ArrayList<Double> distanceList){
		final int SAMPLES = 4;
		double[] lastSAMPLESDifference = new double[SAMPLES];
		double[] lastSAMPLES = new double[SAMPLES];
		double average = 0;
		double min = 99999999;
		for(int i = SAMPLES-1; i>=0; --i){
			int indexBigger = distanceList.size()-i-1;
			int indexSmaller = indexBigger-1;
			lastSAMPLESDifference[i] = distanceList.get(indexBigger) - distanceList.get(indexSmaller);
			lastSAMPLES[i] =  distanceList.get(i);
			min = Math.min(min,lastSAMPLESDifference[i]);
		}

		return true;
	}

	private int mapSeqIDtoMinimumNumberOfSamplesTotake(double seqID){
		return 40;
	}

	public static ArrayList<String> loadFolds(File file, int line)
	{
		ArrayList<String> list = new ArrayList<String>();
		try
		{
			BufferedReader buffer = new BufferedReader(new FileReader(file));
			String textline = null;
			int lines = 0;
			while((textline = buffer.readLine()) != null)
			{
				if((lines - line) % 5 == 0)
				{
					list.add(textline);
				}

				lines++;
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
		return list;
	}
}
