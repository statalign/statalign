package statalign.postprocess.plugins;

import java.awt.BorderLayout;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
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
import com.ppfold.main.PPfoldMain;

import statalign.base.InputData;
import statalign.base.State;
import statalign.base.Utils;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.AlignmentGUI;
import statalign.postprocess.utils.Mapping;
import statalign.postprocess.utils.RNAFoldingTools;

public class PPFold  extends statalign.postprocess.Postprocess  {
	//variables from ppfoldmain

	static private Progress progress= NullProgress.INSTANCE;; //Progressbar; either NullActivity (if no GUI), or the PPfoldProgressBar (if GUI) 

	//end of variables from ppfoldmain


	public String title;	
	public int frequency = 5;
	JPanel pan = new JPanel(new BorderLayout());
	AlignmentGUI gui;
	//private boolean sampling = true;

	CurrentAlignment curAlig;

	ColumnNetwork network;
	Column firstVector, lastVector;
	int sizeOfAlignments;

	int[] firstDescriptor; 
	String t[][];
	String[] sequences;
	String[] viterbialignment;
	int d;
	String refSeq;



	float[][] summedArray;
	int noSambles;

	public PPFold(){
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
		refSeq = input.seqs.sequences.get(0).replaceAll("-","");
		d = refSeq.length();
		
		
		pan.removeAll();
		title = input.title;
		JScrollPane scroll = new JScrollPane();
		scroll.setViewportView(gui = new AlignmentGUI(title,input.model));//, mcmc.tree.printedAlignment()));
		pan.add(scroll, BorderLayout.CENTER);
		pan.getParent().validate();
		
		sizeOfAlignments = (mcmc.tree.vertex.length+1)/2;
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
		for(int i = 0; i < t.length; i++){
			t[i] = curAlig.leafAlignment[i].split("\t");
		}
		Arrays.sort(t, compStringArr);

		int[] previousDescriptor = firstDescriptor;

		int i, j, len = t[0][1].length();
		for(j = 0; j < len; j++){
			int[] nextDescriptor =  new int[sizeOfAlignments];
			boolean allGap = true;
			for(int k = 0; k < sizeOfAlignments; k++){
				if(t[k][1].charAt(j) == '-')
					nextDescriptor[k] = ColumnKey.colNext(previousDescriptor[k]);
				else {
					nextDescriptor[k] = ColumnKey.colNext(previousDescriptor[k])+1;
					allGap = false;
				}
			}
			if(!allGap)
				network.add(nextDescriptor);//[j]);

			previousDescriptor = nextDescriptor;
		}//j (length of alignments)

		if(no == 0) {		// add last vector once only
			int[] lastDescriptor =  new int[sizeOfAlignments];
			for(j = 0; j < sizeOfAlignments; j++){
				lastDescriptor[j] = ColumnKey.colNext(previousDescriptor[j])+1;
			}
			lastVector = network.add(lastDescriptor);
		}
		if(no == 0 || (total-1-no) % frequency == 0) {
			network.updateViterbi(no+1);
			//System.out.println("sequences first: "+sequences);

			if(sequences == null) {
				sequences = new String[sizeOfAlignments];
				for(i = 0; i < sizeOfAlignments; i++){
					sequences[i] = "";
					for(j = 0; j < len; j++){
						if(t[i][1].charAt(j) != '-'){
							sequences[i] += t[i][1].charAt(j);
						}
					}
				}
			}
			for(i = 0; i < sizeOfAlignments; i++)
				viterbialignment[i] = "";
			Column actualVector = lastVector.viterbi;
			ArrayList<Integer> posteriorList = new ArrayList<Integer>();
			while(!actualVector.equals(firstVector)){
				int[] desc = actualVector.key.desc;
				posteriorList.add(new Integer(actualVector.count));
				for(i = 0; i < desc.length; i++){
					if((desc[i] & 1) == 0){
						viterbialignment[i] = "-"+viterbialignment[i];
					}
					else{
						viterbialignment[i] = sequences[i].charAt(desc[i] >> 1) + viterbialignment[i];
					}
				}
				actualVector = actualVector.viterbi;

			}

			if (no==0){
				summedArray = new float[d][d];
				for(i = 0; i<d; ++i){
					for(j = 0; j<d; ++j){
						summedArray[i][j] = 0;
					}
				}
			}


			try {
				List<String> lines = new ArrayList<String>();			
				for (i = 0; i<sequences.length; ++i){
					System.out.println(">" + t[i][0]);
					lines.add(">" + t[i][0].trim());
					System.out.println(t[i][1]);
					lines.add(t[i][1]);
				}	
				
				
				Alignment align = AlignmentReader.readAlignmentFromStringList(lines);

				Tree tree = null;


				BufferedReader paramFileReader = null;
				Parameters param; 
				File file = new File("res/matrices.in");
				paramFileReader = new BufferedReader(new InputStreamReader(	new FileInputStream(file)));

				param = Parameters.readParam(paramFileReader);


				List<ExtraData> extradata = new ArrayList<ExtraData>();
				float[][] fun = PPfoldMain.fold(progress, align.getSequences(), align.getNames(), tree, param, extradata);
				
				float[][] projectFun = Mapping.projectMatrix(t[0][1], fun, '-');
				
				for(i = 0; i<d; ++i){
					for(j = 0; j<d; ++j){
						summedArray[i][j] += projectFun[i][j];
					}
				}
				
				noSambles += 1;
				PPfoldMain.setfoldingfinished(true);

				




				/*
					for(i = 0; i<matrix.length; ++i){
						for(j = 0; j<matrix.length; ++j){
						System.out.print(matrix[i][j]);
						}
						System.out.println();
					}
				 */

			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	

	@Override
	public void afterLastSample(){
		double[][] doubleSummedArray  = new double[d][d];
		for(int i = 0; i<d; ++i){
			for(int j = 0; j<d; ++j){
				doubleSummedArray[i][j] = (double) summedArray[i][j]/noSambles;
			}
		}
	
	
	RNAFoldingTools rnaTools = new RNAFoldingTools();
	
	
	int[] finalmatrix = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(doubleSummedArray);
	String yeah = RNAFoldingTools.getDotBracketStringFromPairedSites(finalmatrix);
	
	System.out.println(finalmatrix);
	System.out.println(yeah);
	
	}





}

