package statalign.postprocess.plugins;

import java.awt.BorderLayout;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import com.ppfold.algo.AsynchronousJobExecutor;
import com.ppfold.algo.AsynchronousJobExecutorThreadPool;
import com.ppfold.algo.ExportTools;
import com.ppfold.algo.FoldingProject;
import com.ppfold.algo.MatrixTools;
import com.ppfold.algo.NeighbourJoining;
import com.ppfold.algo.NullProgress;
import com.ppfold.algo.Parameters;
import com.ppfold.algo.Progress;
import com.ppfold.algo.ResultBundle;
import com.ppfold.algo.Tree;
import com.ppfold.algo.extradata.ExtraData;
import com.ppfold.main.Alignment;
import com.ppfold.main.AlignmentReader;
import com.ppfold.main.DataInfo;
import com.ppfold.main.NewickReader;
import com.ppfold.main.PPfoldMain;

import statalign.base.InputData;
import statalign.base.Mcmc;
import statalign.base.State;
import statalign.base.Utils;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.AlignmentGUI;
import statalign.postprocess.gui.TestGUI;
import statalign.postprocess.utils.Mapping;
import statalign.postprocess.utils.RNAFoldingTools;

public class PPFold extends statalign.postprocess.Postprocess {
	// variables from ppfoldmain

	static private Progress progress = NullProgress.INSTANCE;; // Progressbar;
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
	TestGUI gui;
	// private boolean sampling = true;

	CurrentAlignment curAlig;

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

	float[][] summedArray;
	float[] summedSingleBaseProb;
	float[][] probMatrix;
	
	int noSambles;

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
		return "PPFold Plugin";
	}

	@Override
	public void setSampling(boolean enabled) {
		sampling = enabled;
	}

	@Override
	public String[] getDependences() {
		return new String[] { "statalign.postprocess.plugins.CurrentAlignment" };
	}

	@Override
	public void refToDependences(Postprocess[] plugins) {
		curAlig = (CurrentAlignment) plugins[0];
	}

	static Comparator<String[]> compStringArr = new Comparator<String[]>() {
		public int compare(String[] a1, String[] a2) {
			return a1[0].compareTo(a2[0]);
		}
	};

	@Override
	public void beforeFirstSample(InputData input) {
		int maxLength = 0;
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

		pan.removeAll();
		title = input.title;
		JScrollPane scroll = new JScrollPane();
		scroll.setViewportView(gui = new TestGUI(this));// ,
																			// mcmc.tree.printedAlignment()));
		pan.add(scroll, BorderLayout.CENTER);
		pan.getParent().validate();
		
		gui.changeDimension(d);
		sizeOfAlignments = (mcmc.tree.vertex.length + 1) / 2;
		noSambles = 0;

		t = new String[sizeOfAlignments][];
		sequences = null;
		viterbialignment = new String[sizeOfAlignments];

		network = new ColumnNetwork();

		firstDescriptor = new int[sizeOfAlignments];
		Arrays.fill(firstDescriptor, -1);
		firstVector = network.add(firstDescriptor);
		lastVector = null;

		viterbialignment = new String[sizeOfAlignments];
	}

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
				summedArray = new float[d][d];
				for (i = 0; i < d; ++i) {
					for (j = 0; j < d; ++j) {
						summedArray[i][j] = 0;
					}
				}
				summedSingleBaseProb = new float[d];
			}

			try {
				List<String> lines = new ArrayList<String>();
				for (i = 0; i < sequences.length; ++i) {
					System.out.println(">" + t[i][0]);
					lines.add(">" + t[i][0].trim());
					System.out.println(t[i][1]);
					lines.add(t[i][1]);
				}

				Alignment align = AlignmentReader
						.readAlignmentFromStringList(lines);

				// Tree tree = null;

				Tree tree = getPPfoldTree(mcmc);

				BufferedReader paramFileReader = null;
				Parameters param;
				File file = new File("res/matrices.in");
				paramFileReader = new BufferedReader(new InputStreamReader(
						new FileInputStream(file)));

				param = Parameters.readParam(paramFileReader);

				List<ExtraData> extradata = new ArrayList<ExtraData>();
				float[][] fun = PPfoldMain.fold(progress, align.getSequences(),
						align.getNames(), tree, param, extradata);
				float[] singleBaseProb = new float[fun.length];
				for (int x = 0; x < fun.length; x++) {
					singleBaseProb[x] = 1;
					for (int y = 0; y < fun[0].length; y++) {
						singleBaseProb[x] -= fun[x][y];
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
						PPFold.getSequenceByName(t, refSeqName), fun, '-');
				// float [] projectSingleBaseProb = Mapping.projectArray(mapSeq,
				// singleBaseProb, '-');
				float[] projectSingleBaseProb = Mapping.projectarray(
						PPFold.getSequenceByName(t, refSeqName),
						singleBaseProb, '-');
				
				// normalise project matrix
				for(int x = 0 ; x < projectFun.length ; x++)
				{
					double rowMatrSum = 0;
					for(int y = 0 ; y < projectFun[0].length ; y++)
					{
						rowMatrSum += projectFun[x][y];
					}
					
					double factor = rowMatrSum + projectSingleBaseProb[x];
					System.out.println("F:"+factor);
					for(int y = 0 ; y < projectFun[0].length ; y++)
					{
						projectFun[x][y] = (float)(projectFun[x][y] / factor);
						projectSingleBaseProb[x] /= factor;
					}
				}

				/*
				System.out.println("D=" + d);
				System.out.println(summedArray.length);
				System.out.println(projectFun.length);*/
				probMatrix = new float[d][d];
				for (i = 0; i < d; ++i) {
					summedSingleBaseProb[i] += projectSingleBaseProb[i];
					for (j = 0; j < d; ++j) {
						summedArray[i][j] += projectFun[i][j];
						probMatrix[i][j] = summedArray[i][j]/(float)(noSambles+1);
						gui.setMatrix(probMatrix);
						gui.repaint();
					}
				}

				noSambles += 1;
				
			
				RNAFoldingTools rnaTools = new RNAFoldingTools();
				String seq = this.getSequenceByName(t, this.refSeqName).replaceAll("-", "");
				int[] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(probMatrix);
				System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites));
				
				Structure.updateMatrix(probMatrix);
				
				PPfoldMain.setfoldingfinished(true);
				
				
				/*float[][] probMatrix = new float[d][d];
				for(int x = 0; x < d; x++) {
					for(int y = 0; y < d; y++) {
						probMatrix[x][y] = summedArray[x][y]/noSambles;
					}
				}*/
				//if(sampling) {
				gui.clear();
				gui.changeDimension(d*TestGUI.OFFSET);
				gui.setMatrix(probMatrix);
				gui.repaint();
				//}

				/*
				 * for(i = 0; i<matrix.length; ++i){ for(j = 0; j<matrix.length;
				 * ++j){ System.out.print(matrix[i][j]); } System.out.println();
				 * }
				 */

			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	@Override
	public void afterLastSample() {
		double[][] doubleSummedArray = new double[d][d];
		double[] doubleSingleBaseProb = new double[d];
		for (int i = 0; i < d; ++i) {
			doubleSingleBaseProb[i] = summedSingleBaseProb[i]
					/ (double) noSambles;
			for (int j = 0; j < d; ++j) {
				doubleSummedArray[i][j] = (double) summedArray[i][j]
						/ (double) noSambles;
			}
			// System.out.println();
		}

		RNAFoldingTools rnaTools = new RNAFoldingTools();

		int[] pairedSites = rnaTools
				.getPosteriorDecodingConsensusStructureMultiThreaded(doubleSummedArray);
		// int[] finalmatrix =
		// rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(doubleSummedArray);
		
		RNAFoldingTools.writeMatrix(RNAFoldingTools.getDoubleMatrix(probMatrix), new File("prob.matrix"));
		
		
		RNAFoldingTools.writeMatrix(doubleSummedArray, new File("bp.matrix"));
		System.out.println("num samples" + noSambles);
		//RNAFoldingTools.writeMatrix(summedArray, new File("bp2.matrix"));
		double[] singleBaseProb = RNAFoldingTools
				.getSingleBaseProb(doubleSummedArray);
		saveResult(refSeqGapped, pairedSites, doubleSummedArray,
				singleBaseProb, new File("result.txt"));
	}

	public static String getSequenceByName(String[][] sequences, String name) {
		for (int i = 0; i < sequences.length; i++) {
			if (sequences[i][0].equals(name)) {
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
	
}
