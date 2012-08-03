package statalign.postprocess.plugins;

import java.awt.BorderLayout;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
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

import statalign.base.InputData;
import statalign.base.Mcmc;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.base.Utils;
import statalign.distance.Distance;
import statalign.postprocess.Postprocess;
import statalign.postprocess.plugins.benchmarks.Benchmarks;
import statalign.postprocess.gui.PPFoldGUI;
import statalign.postprocess.utils.Mapping;
import statalign.postprocess.utils.RNAFoldingTools;

public class PPFold extends statalign.postprocess.Postprocess {

	public static void main(String [] args)
	{
		int [] pairedSites = RNAFoldingTools.getPosteriorDecodingConsensusStructure(RNAFoldingTools.loadMatrix(new File(System.getProperty("user.home") + "/Dropbox/RNA and StatAlign/static/bp.matrix")));
		System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites));
	}

	// variables from ppfoldmain

	static private Progress progress = NullProgress.INSTANCE; // Progressbar;
	// either
	// NullActivity
	// (if no GUI),
	// or the
	// PPfoldProgressBar
	// (if GUI)	

	// end of variables from ppfoldmain

	public String title;
	//public int frequency = 5;
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
	String refSeqGapped;

	double[][] summedBasePairProbMatrix;
	double[] summedSingleBaseProb;
	float[][] probMatrix;

	int noSamples;

	double [][] weightedBasePairProb;
	double beta = 10;
	double weightedSum = 0;
	double firstLikelihood = 0;
	double posteriorProbabilityAvg = 0;
	ArrayList<String> tempAlignment;
	ArrayList<Double> logLikelihood;

	RNAFoldingTools rnaTools = new RNAFoldingTools();

	String outDir = "output/";

	public PPFold() {
		screenable = true;
		outputable = true;
		postprocessable = true;
		postprocessWrite = true;
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
		public int compare(String[] a1, String[] a2) {
			return a1[0].compareTo(a2[0]);
		}
	};



	List<ExtraData> extradata = new ArrayList<ExtraData>();
	Parameters param = null;

	@Override
	public void beforeFirstSample(InputData input) {
		int maxLength = 0;
		title = input.title;
		// refSeq = input.seqs.sequences.get(0);
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

		// System.out.println("using seq no. " + seqNo);

		d = refSeq.length();


		if(show) {
			pan.removeAll();
			title = input.title;
			JScrollPane scroll = new JScrollPane();
			scroll.getViewport().add(gui = new PPFoldGUI(title, scroll));

			if(gui.getPreferredSize().getHeight() > gui.getHeight()) {
				scroll.createVerticalScrollBar();
			}

			if(gui.getPreferredSize().getWidth() > gui.getWidth()) {
				scroll.createHorizontalScrollBar();
			}

			pan.add(scroll, BorderLayout.CENTER);
			pan.getParent().validate();

		}

		if(gui != null)
		{
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
	}

	public double saveMPDToFile(String [] fastaAlignment, File outFile)
	{
		double ppfoldReliability = 0;
		List<String> lines = new ArrayList<String>();
		for (int i = 0; i < fastaAlignment.length; i++) {			
			lines.add(fastaAlignment[i]);
		}

		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outDir+title+"_mpd.fas"));
			for (int i = 0; i < fastaAlignment.length; i++) {
				buffer.write(fastaAlignment[i]+"\n");
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}


		try {

			Alignment align = AlignmentReader.readAlignmentFromStringList(lines);
			List<String> sequences = align.getSequences();
			List<String> seqNames = align.getNames();
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

			ResultBundle mpdResult = PPfoldMain.fold(progress, align.getSequences(), align.getNames(), null, param, extradata);
			RNAFoldingTools.writeToFile(new File(outDir + title+"_entropy_mpd.txt"), "mpd\t"+mpdResult.entropyVal+"\t"+mpdResult.entropyPercOfMax + "\t"+mpdResult.entropyMax, true);
			float [][] basePairProb = mpdResult.finalmatrix;			
			System.out.println("MPD entropy" + mpdResult.entropyVal+"\t"+mpdResult.entropyPercOfMax+"\t"+mpdResult.entropyMax);

			int [] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(basePairProb);
			double ppfoldReliability2 = RNAFoldingTools.calculatePPfoldReliabilityScore(pairedSites, RNAFoldingTools.getDoubleMatrix(basePairProb));
			ppfoldReliability = mpdResult.getPPfoldReliability();
			System.out.println("Reliabilities: " + ppfoldReliability+"\t"+ppfoldReliability2);
			int [] projectedPairedSites = Benchmarks.projectPairedSites(refSeqGapped, pairedSites);
			System.out.println("PPfold reliability score (MPD) = " + ppfoldReliability);
			System.out.println("Improved reliability score 1 (MPD) = " + (posteriorProbabilityAvg*ppfoldReliability));
			double proxySimilarity = getProxySimilarityFromAvgPosterior(posteriorProbabilityAvg);
			System.out.println("Improved reliability score 2 (MPD) = " + (proxySimilarity*ppfoldReliability));


			try
			{
				BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile));
				buffer.write(refSeqGapped+"\n");
				buffer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(projectedPairedSites)+"\n");
				buffer.close();
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		return ppfoldReliability;
	}

	ArrayList<AlignmentData> alignments = new ArrayList<AlignmentData>();
	String p2 = "";
	@Override
	public void newSample(State state, int no, int total) {		

		for (int i = 0; i < t.length; i++) {
			t[i] = curAlig.leafAlignment[i].split("\t");
		}
		Arrays.sort(t, compStringArr);
	
		int[] previousDescriptor = firstDescriptor;

		int i, j, len = t[0][1].length();
		for (j = 0; j < len; j++) {
			int[] nextDescriptor = new int[sizeOfAlignments];
			boolean allGap = true;
			for (int k = 0; k < sizeOfAlignments; k++) {
				if (t[k][1].charAt(j) == '-')
					nextDescriptor[k] = ColumnKey
					.colNext(previousDescriptor[k]);
				else {
					nextDescriptor[k] = ColumnKey
							.colNext(previousDescriptor[k]) + 1;
					allGap = false;
				}
			}
			if (!allGap)
				network.add(nextDescriptor);// [j]);

			previousDescriptor = nextDescriptor;
		}// j (length of alignments)

		if (no == 0) { // add last vector once only
			int[] lastDescriptor = new int[sizeOfAlignments];
			for (j = 0; j < sizeOfAlignments; j++) {
				lastDescriptor[j] = ColumnKey.colNext(previousDescriptor[j]) + 1;
			}
			lastVector = network.add(lastDescriptor);

		}
		if (no == 0 || 1 < 2) {
			network.updateViterbi(no + 1);
			// System.out.println("sequences first: "+sequences);

			if (sequences == null) {
				sequences = new String[sizeOfAlignments];
				for (i = 0; i < sizeOfAlignments; i++) {
					sequences[i] = "";
					for (j = 0; j < len; j++) {
						if (t[i][1].charAt(j) != '-') {
							sequences[i] += t[i][1].charAt(j);
						}
					}
				}
			}
			for (i = 0; i < sizeOfAlignments; i++)
				viterbialignment[i] = "";
			Column actualVector = lastVector.viterbi;
			ArrayList<Integer> posteriorList = new ArrayList<Integer>();
			while (!actualVector.equals(firstVector)) {
				int[] desc = actualVector.key.desc;
				posteriorList.add(new Integer(actualVector.count));
				for (i = 0; i < desc.length; i++) {
					if ((desc[i] & 1) == 0) {
						viterbialignment[i] = "-" + viterbialignment[i];
					} else {
						viterbialignment[i] = sequences[i].charAt(desc[i] >> 1)
								+ viterbialignment[i];
					}
				}
				actualVector = actualVector.viterbi;

			}

			if (no == 0) {
				summedBasePairProbMatrix = new double[d][d];
				weightedBasePairProb = new double[d][d];
				for (i = 0; i < d; ++i) {
					for (j = 0; j < d; ++j) {
						summedBasePairProbMatrix[i][j] = 0;
						weightedBasePairProb[i][j] = 0;
					}
				}
				summedSingleBaseProb = new double[d];
			}

			try {
				List<String> lines = new ArrayList<String>();
				for (i = 0; i < sequences.length; ++i) {
					lines.add(">" + t[i][0].trim());
					lines.add(t[i][1]);
				}
				
				Alignment align = AlignmentReader.readAlignmentFromStringList(lines);

				Tree tree = getPPfoldTree(mcmc);
				
				//writes the distance of the first sampleq to a file
				BufferedWriter buffer2 = new BufferedWriter(new FileWriter(outDir+title+"_sample1SeqID"));
				System.out.println(Distance.sequenceSimilarityScore(align.getSequences()));
				buffer2.write(""+Distance.sequenceSimilarityScore(align.getSequences()));
				buffer2.close();
				
				
				//AlignmentData d = new AlignmentData();
				ResultBundle sampleResult = PPfoldMain.fold(progress, align.getSequences(),	align.getNames(), tree, param, extradata);
				System.out.println("Sample "+no+ " entropy" + sampleResult.entropyVal+"\t"+sampleResult.entropyPercOfMax+"\t"+sampleResult.entropyMax);
				RNAFoldingTools.writeToFile(new File(outDir + title+"_entropy_samples.txt"), no+"\t"+sampleResult.entropyVal+"\t"+sampleResult.entropyPercOfMax + "\t"+sampleResult.entropyMax, true);
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
					// System.out.println(k+"\t"+t[k][1]);
				}

				// System.out.println(t[seqNo][1]);	
				// System.out.println(refSeq);
				// String mapSeq =
				// RNAFoldingTools.getReferenceSequence(sequences,
				// refSeq.length());
				// System.out.println("REFSEQ: " + mapSeq);
				// float[][] projectFun = Mapping.projectMatrix(mapSeq, fun,
				// '-');
				float[][] projectFun = Mapping.projectMatrix(
						PPFold.getSequenceByName(t, refSeqName), basePairProb, '-');
				// float [] projectSingleBaseProb = Mapping.projectArray(mapSeq,
				// singleBaseProb, '-');
				float[] projectSingleBaseProb = Mapping.projectarray(
						PPFold.getSequenceByName(t, refSeqName),
						singleBaseProb, '-');

				// normalise projected matrix
				for(int x = 0 ; x < projectFun.length ; x++)
				{
					double rowMatrSum = 0;
					for(int y = 0 ; y < projectFun[0].length ; y++)
					{
						rowMatrSum += projectFun[x][y];
					}
					double factor = rowMatrSum + projectSingleBaseProb[x];
					factor = 1;
					//
					System.out.println("F:"+factor+"\t"+rowMatrSum+"\t"+projectSingleBaseProb[x]);
					for(int y = 0 ; y < projectFun[0].length ; y++)
					{
						projectFun[x][y] = (float)(projectFun[x][y] / factor);
					}
					projectSingleBaseProb[x] /= factor;
				}

				double alignmentLogLikelihood = mcmc.mcmcStep.newLogLike;
				if(noSamples == 0)
				{
					firstLikelihood = mcmc.mcmcStep.newLogLike;
					weightedSum  = 0;
					try {
						BufferedWriter bufferClear = new BufferedWriter(new FileWriter(outDir + title + "likelihoods.txt", false));
						bufferClear.write("");
						bufferClear.close();
					}catch (IOException ex){
						ex.printStackTrace();
					}

				}
				try
				{

					BufferedWriter buffer = new BufferedWriter(new FileWriter(outDir + title + "likelihoods.txt", true));
					buffer.write(noSamples+"\t"+(alignmentLogLikelihood - firstLikelihood)+"\n");

					buffer.close();

					//System.out.println(noSambles+"\t"+mcmc.mcmcStep.newLogLike);
				}
				catch(IOException ex)
				{
					ex.printStackTrace();
				}

				RNAFoldingTools rnaTools = new RNAFoldingTools();
				boolean append = true;
				if(noSamples == 0)
				{
					append = false;
				}
				PPFold.appendFolds(new File(outDir+title+".folds"), noSamples+"", PPFold.getSequenceByName(t, refSeqName),rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(basePairProb), rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(projectFun), append);





				/*
				System.out.println("D=" + d);
				System.out.println(summedArray.length);
				System.out.println(projectFun.length);*/
				probMatrix = new float[d][d];
				double weight = Math.pow(firstLikelihood / alignmentLogLikelihood, beta);

				for (i = 0; i < d; ++i) {
					summedSingleBaseProb[i] += projectSingleBaseProb[i];
					for (j = 0; j < d; ++j) {
						summedBasePairProbMatrix[i][j] += projectFun[i][j];

						probMatrix[i][j] = (float)(summedBasePairProbMatrix[i][j]/(double)(noSamples+1));
						if(gui != null)
						{
							gui.setMatrix(probMatrix);
							gui.repaint();
						}

						weightedBasePairProb[i][j] += projectFun[i][j]*weight;

					}
				}
				weightedSum += weight;

				BufferedWriter buffer = new BufferedWriter(new FileWriter("weights.txt", true));
				buffer.write(noSamples+"\t"+weightedSum+"\t"+weight+"\n");
				buffer.close();

				//System.out.println(noSambles+"\t"+mcmc.mcmcStep.newLogLike);
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}

			noSamples += 1;

			//RNAFoldingTools rnaTools = new RNAFoldingTools();
			//String seq = this.getSequenceByName(t, this.refSeqName).replaceAll("-", "");
			//RNAFoldingTools rnaTools = new RNAFoldingTools();
			int[] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(probMatrix);
			p2 = RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites);
			System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites));

			Structure.updateMatrix(probMatrix);

			PPfoldMain.setfoldingfinished(true);

			if(gui != null)
			{
				gui.changeDimension(d*PPFoldGUI.OFFSET);
				gui.setMatrix(probMatrix);
				gui.repaint();
			}

		}

		// save alignment sample information
		if(no == 0)
		{
			String [] [] inputAlignment = new String[mpdAlignment.input.seqs.sequences.size()][2];

			for(int k = 0 ; k < inputAlignment.length ; k++)
			{
				inputAlignment[k][0] = mpdAlignment.input.seqs.seqNames.get(k);
				inputAlignment[k][1] = mpdAlignment.input.seqs.sequences.get(k);
			}
			Arrays.sort(inputAlignment, compStringArr);

			/*//AlignmentData al =  new AlignmentData();
			for(int k = 0 ; k < t.length ; k++)
			{
				al.sequences.add(t[k][1]);
				al.names.add(t[k][0]);
			}
			alignments.add(al);*/

			appendAlignment("reference", inputAlignment, new File(outDir+mpdAlignment.input.title+".samples"), false);
			appendAlignment(no+"", t, new File(outDir+mpdAlignment.input.title+".samples"), true);
			RNAFoldingTools.writeToFile(new File(outDir +title+"_Distance.txt"),"BurnIn "+new Integer(mcmc.mcmcpars.burnIn).toString() + "\t" +
					"cycles "+new Integer(mcmc.mcmcpars.cycles).toString() + "\t" +
					"sampRate "+new Integer(mcmc.mcmcpars.sampRate).toString() + "\t" + 
					"Filename "+ mpdAlignment.input.title + "\t", false);


			tempAlignment = new ArrayList<String>();
			for(int k = 0 ; k < t.length ; k++)
			{
				tempAlignment.add(t[k][1]);
			}
			//al.


			RNAFoldingTools.writeToFile(new File(outDir +title+"_entropy_fuzzy_obs.txt"),"", false);
			RNAFoldingTools.writeToFile(new File(outDir +title+"_entropy_fuzzy_exp.txt"),"", false);
			RNAFoldingTools.writeToFile(new File(outDir +title+"_entropy_samples.txt"),"", false);
			RNAFoldingTools.writeToFile(new File(outDir +title+"_entropy_mpd.txt"),"", false);
			alignments = new ArrayList<AlignmentData>();

			//for the MCMC part (part 4)



		}
		else
		{
			appendAlignment(no+"", t, new File(outDir+mpdAlignment.input.title+".samples"), true);
			//alignments.clear();
			//Arrays.sort(t, compStringArr);
			AlignmentData al =  new AlignmentData();
			for(int k = 0 ; k < t.length ; k++)
			{
				al.sequences.add(t[k][1]);
				al.names.add(t[k][0]);
			}
			alignments.add(al);


			FuzzyAlignment fuzzyAlignment = FuzzyAlignment.getFuzzyAlignmentAndProject(alignments, refSeqName);

			double dist = Distance.AMA(tempAlignment, al.sequences);

			/*if(dist < CONST_DIST){
				tempAlignment.clear();
				tempAlignment.addAll(al.sequences);
				RNAFoldingTools.writeToFile(new File(outDir +title+"_Distance.txt"),new Double(dist).toString(), true);
				RNAFoldingTools.writeToFile(new File(outDir +title+"_Distance.txt"),"1", true);
			}
			else{*/
			RNAFoldingTools.writeToFile(new File(outDir +title+"_Distance.txt"),new Double(dist).toString(), true);
			//}

			List<ExtraData> extradata = new ArrayList<ExtraData>();
			try {
				System.out.println(fuzzyAlignment.toString());
				Tree tree = null;
				//Tree tree = getPPfoldTree(mcmc);
				System.out.println("Tree:");
				getPPfoldTree(mcmc).print();
				ResultBundle fuzzyResultExp = PPfoldMain.foldFuzzyAlignment(progress, fuzzyAlignment, tree, param, extradata, true);
				ResultBundle fuzzyResultObs = PPfoldMain.foldFuzzyAlignment(progress, fuzzyAlignment, tree, param, extradata, false);
				int [] fuzzyPairedSitesExp = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(fuzzyResultExp.finalmatrix);
				int [] fuzzyPairedSitesObs = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(fuzzyResultObs.finalmatrix);

				System.out.println("R0:"+p2);
				System.out.println(refSeq);				
				System.out.println("FZ:"+RNAFoldingTools.getDotBracketStringFromPairedSites(fuzzyPairedSitesExp));
				System.out.println("Fuzzy entropy: " + fuzzyResultExp.entropyVal+"\t"+fuzzyResultExp.entropyPercOfMax+"\t"+fuzzyResultExp.entropyMax);
				RNAFoldingTools.writeToFile(new File(outDir + title+"_entropy_fuzzy_exp.txt"), no+"\t"+fuzzyResultExp.entropyVal+"\t"+fuzzyResultExp.entropyPercOfMax + "\t"+fuzzyResultExp.entropyMax, true);
				RNAFoldingTools.writeToFile(new File(outDir + title+"_entropy_fuzzy_obs.txt"), no+"\t"+fuzzyResultObs.entropyVal+"\t"+fuzzyResultObs.entropyPercOfMax + "\t"+fuzzyResultObs.entropyMax, true);
				PPFold.appendFolds(new File(outDir+title+".folds_e_exp"), "exp", PPFold.getSequenceByName(t, refSeqName),fuzzyPairedSitesExp,fuzzyPairedSitesExp, false);
				PPFold.appendFolds(new File(outDir+title+".folds_e_obs"), "obs", PPFold.getSequenceByName(t, refSeqName),fuzzyPairedSitesExp,fuzzyPairedSitesObs, false);
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	int stepCounter = 0;


	@Override
	public void afterLastSample() {
		double[][] doubleSummedArray = new double[d][d];
		double[] doubleSingleBaseProb = new double[d];
		for (int i = 0; i < d; ++i) {
			doubleSingleBaseProb[i] = summedSingleBaseProb[i]
					/ (double) noSamples;
			for (int j = 0; j < d; ++j) {
				doubleSummedArray[i][j] = (double) summedBasePairProbMatrix[i][j]
						/ (double) noSamples;
			}
		}

		for (int i = 0; i < d; ++i) {
			double sum = doubleSingleBaseProb[i];
			for (int j = 0; j < d; ++j) {
				sum += doubleSummedArray[i][j];
			}
			System.out.println(i+" > " + sum);
		}

		RNAFoldingTools rnaTools = new RNAFoldingTools();

		int[] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(doubleSummedArray);
		double statalignPpfoldReliablityScore = RNAFoldingTools.calculatePPfoldReliabilityScore(pairedSites, doubleSummedArray);
		//RNAFoldingTools.writeMatrix(RNAFoldingTools.getDoubleMatrix(probMatrix), new File("prob.matrix"));


		RNAFoldingTools.writeMatrix(doubleSummedArray, new File("bp.matrix"));
		System.out.println("num samples" + noSamples);
		//RNAFoldingTools.writeMatrix(summedArray, new File("bp2.matrix"));
		double[] singleBaseProb = RNAFoldingTools.getSingleBaseProb(doubleSummedArray);
		saveResult(refSeqGapped, pairedSites, doubleSummedArray, singleBaseProb, new File(title+".dat.res"));

		for (int i = 0; i < d; ++i) {
			for (int j = 0; j < d; ++j) {
				this.weightedBasePairProb[i][j] /= weightedSum;
			}
		}
		RNAFoldingTools.writeMatrix(weightedBasePairProb, new File(outDir + "bp_log.matrix"));
		saveResult(refSeqGapped, rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(weightedBasePairProb), weightedBasePairProb, RNAFoldingTools.getSingleBaseProb(weightedBasePairProb), new File(outDir +title+".dat.res.weighted"));


		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(new File(outDir+mpdAlignment.input.title+".samples"), true));
			buffer.write("%posteriors\n");
			for(int i = 0 ; i < mpdAlignment.decoding.length ; i++)
			{
				buffer.write(mpdAlignment.decoding[i]+"\n");
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
		posteriorProbabilityAvg = 0;
		for(int i = 0 ; i < mpdAlignment.decoding.length ; i++)
		{
			posteriorProbabilityAvg += mpdAlignment.decoding[i];
		}
		posteriorProbabilityAvg /= (double)(mpdAlignment.decoding.length);
		System.out.println("Posterior probability avg: " + posteriorProbabilityAvg);
		System.out.println("Structure reliability score = " + statalignPpfoldReliablityScore);

		double ppfoldReliablityScore = RNAFoldingTools.calculatePPfoldReliabilityScore(pairedSites, doubleSummedArray);
		System.out.println("Structure reliability score = " + ppfoldReliablityScore);
		double proxySim = getProxySimilarityFromAvgPosterior(posteriorProbabilityAvg);
		System.out.println("Proxy similiarity = " + proxySim);
		ArrayList<String> sequences = new ArrayList<String>();
		for(int i = 0 ; i < mpdAlignment.alignment.length ; i++)
		{
			sequences.add(mpdAlignment.alignment[i]);
		}			
		System.out.println("Proxy distance = " + Distance.amaScoreToMultiDistance(sequences,proxySim));
		double improvedReliabilityScore1 = posteriorProbabilityAvg*statalignPpfoldReliablityScore;		
		System.out.println("Improved reliability score 1 = " + improvedReliabilityScore1);
		double improvedReliabilityScore2 = proxySim*statalignPpfoldReliablityScore;		
		System.out.println("Improved reliability score 2 = " + improvedReliabilityScore2);

		// save mpd alignment
		appendAlignment("mpd", mpdAlignment.alignment, new File(outDir+mpdAlignment.input.title+".samples"), true, mpdAlignment.input);
		double ppfoldReliabilityScoreMPD = saveMPDToFile(Utils.alignmentTransformation(mpdAlignment.alignment,
				"Fasta", mpdAlignment.input), new File(outDir+mpdAlignment.input.title+".dat.res.mpd"));
		try
		{

			BufferedWriter buffer = new BufferedWriter(new FileWriter(outDir+mpdAlignment.input.title+".samples", true));
			buffer.write("%logLiklihood\n");
			for(int i = 0; i<logLikelihood.size(); ++i){
				buffer.write(logLikelihood.get(i) + "\n");
			}

			buffer.close();

			//System.out.println(noSambles+"\t"+mcmc.mcmcStep.newLogLike);
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}

		this.saveScores(new File(outDir +"scores.txt"), ppfoldReliablityScore, ppfoldReliabilityScoreMPD);
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
			double proxyDistance = Distance.amaScoreToMultiDistance(sequences,proxySimilarity);
			buffer.write(mpdAlignment.input.title+"\t");
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
			System.out.println(sequences[i]+"\t"+name);
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
