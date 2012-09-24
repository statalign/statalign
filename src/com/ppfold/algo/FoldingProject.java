package com.ppfold.algo;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import statalign.postprocess.utils.RNAFoldingTools;

import com.ppfold.algo.extradata.ExtraData;
import com.ppfold.algo.extradata.ForcedConstraints;

/**
 * After reading and parsing input, algorithm execution starts here. (Note:
 * Columns with many gaps are already removed)
 * 
 * @author Z.Sukosd
 */
public class FoldingProject {
	static final double LOG_TWO = Math.log(2d);

	/**
	 * Makes sure that all leaves of the tree can be found, and creates a list
	 * of leaves for easy access.
	 */
	public static boolean checkInput(Tree tree, List<char[]> columns,
			List<String> names) {
		boolean treematchesnames = true;
		boolean columnshaverightlengths = true;
		boolean noillegalcharacters = true;

		int firstsize = columns.get(0).length;
		for (int col = 0; col < columns.size(); col++) {
			if (columns.get(col).length != firstsize) {
				columnshaverightlengths = false;
				System.err
						.println("Columns in alignment don't have same length!");
			}
			for (int row = 0; row < columns.get(col).length; row++) {
				char thischar = columns.get(col)[row];
				thischar = Character.toLowerCase(thischar);
				if (thischar != 'a' && thischar != 'u' && thischar != 't'
						&& thischar != 'g' && thischar != 'c'
						&& thischar != 'r' && thischar != 'y'
						&& thischar != 's' && thischar != 'w'
						&& thischar != 'k' && thischar != 'm'
						&& thischar != 'b' && thischar != 'd'
						&& thischar != 'h' && thischar != 'v'
						&& thischar != 'n'
						&& !MatrixTools.isGap(thischar)) {
					System.err.println("Illegal character! " + thischar);
					noillegalcharacters = false;
				}
			}
		}

		for (int row = 0; row < columns.get(0).length; row++) {
			// finds the node corresponding to the rownumber
			// rownumber = corresponds to sequences.
			Node node = tree.findSlowlyNodeWithName(names.get(row));
			if (node == null) {
				System.err.println("Can't find node with name "
						+ names.get(row) + "!");
				treematchesnames = false;
			}
		}
		return treematchesnames && columnshaverightlengths
				&& noillegalcharacters;
	}

//	public static ResultBundle fold(int phylojobsnr, int scfgjobsnr, Tree tree,
//			List<char[]> columns, List<String> names, Parameters param,
//			AsynchronousJobExecutor executor, boolean verbose, int execnr,
//			List <ExtraData> extradata_list, boolean diffbp, boolean entropycalc) {
//		try {
//			return fold(NullProgress.INSTANCE, phylojobsnr, scfgjobsnr, tree,
//					columns, names, param, executor, verbose, execnr, extradata_list, diffbp, entropycalc);
//		} catch (InterruptedException e) {
//			// Never happens - we use the NullActivity
//			e.printStackTrace();
//			return null;
//		}
//	}

	public static ResultBundle fold(Progress act, int phylojobsnr,
			int scfgjobsnr, Tree tree, List<char[]> columns,
			List<String> names, Parameters param,
			AsynchronousJobExecutor executor, boolean verbose, int execnr,
			List <ExtraData> extradata_list, boolean diffbp, boolean entropycalc)
			throws InterruptedException {
		
		if(columns.size()<2){
			return ResultBundle.tinyBundle();
		}
		
		double scfg_to_phylo_ratio = 0.095 * columns.size()
				/ columns.get(0).length;
		double phylopart = 0.95 / (scfg_to_phylo_ratio + 1);
		if (phylopart < 0 | phylopart > 1) {
			System.err.println("Time estimation resulted "
					+ "in illegal value for phylogenetic calculations: "
					+ phylopart);
			phylopart = 0.475; // in case the estimation results in a
			// "weird number", just set it to 50-50%.
		}
		double scfgpart = 0.95 - phylopart; // max is 95% because processing
		// must happen before and after
		act.setCurrentActivity("Applying evolutionary model");

		
		//long starttime = System.currentTimeMillis();
		double[][] probmatrix = FoldingProject.createPhyloProb(act
				.getChildProgress(phylopart), phylojobsnr, tree, columns,
				names, columns.size(), param, executor, verbose, execnr);

		//If there are extra data, multiply probmatrix with the relevant values
		
		double[][] probmatrix2 = null;
		if(diffbp&&extradata_list!=null){
			probmatrix2= new double[probmatrix.length][probmatrix[0].length];
			for (int a = 0; a<probmatrix2.length;a++){
				for(int b = 0;b<probmatrix2[0].length; b++){
					probmatrix2[a][b] = probmatrix[a][b]; 
				}
			}
		}
		
		if(extradata_list!=null){
			System.out.println("Number of auxiliary data items: " + extradata_list.size());
			ArrayList<ForcedConstraints> fcList = new ArrayList<ForcedConstraints>();
			for(ExtraData data:extradata_list){
				if(data instanceof ForcedConstraints){
					//If forcing constraints, collect them into one and remove from the list
					fcList.add((ForcedConstraints) data);
				}
			}
			for(ExtraData data:fcList){
				extradata_list.remove(data);
			}
			if(fcList.size()>0){
				System.out.println("Processing " + fcList.size () + " hard constraint datasets");
				ForcedConstraints fcMerged = new ForcedConstraints();
				fcMerged = fcMerged.combinedForcedConstraints(probmatrix.length, fcList);
				extradata_list.add(fcMerged);
			}

		}
		
		System.out.println("Processing all auxiliary data");
		if(extradata_list!=null){
			for(ExtraData extradata:extradata_list){
				for(int i = 0; i<probmatrix.length; i++){
					for(int j = 0; j<i; j++){
						probmatrix[i][j] *= extradata.getProbabilityGivenInnerPaired(i,j);
						probmatrix[j][i] = probmatrix[i][j];
					}
					probmatrix[i][i] *= extradata.getProbabilityGivenUnpaired(i);
				}
			}
		}
		
		if(diffbp&&extradata_list!=null){
			//If inner and outer basepairs are distinguished
			for(ExtraData extradata:extradata_list){
				for(int i = 0; i<probmatrix2.length; i++){
					for(int j = 0; j<i; j++){
						probmatrix2[i][j] *= extradata.getProbabilityGivenOuterPaired(i,j);
						probmatrix2[j][i] = probmatrix2[i][j];
					}
					probmatrix2[i][i] *= extradata.getProbabilityGivenUnpaired(i);
				}
			}
			
		}
		
		//System.out.println("Time in phylogenetic part: " + (System.currentTimeMillis()-starttime));
		// do not allow basepairs less than 4 nucleotides apart.
		for (int i = 0; i < probmatrix.length; i++) {
			for (int j = 1; j < 4; j++) {
				if (i + j < probmatrix.length) {
					probmatrix[i][i + j] = 0;
				}
				if (i - j > 0) {
					probmatrix[i][i - j] = 0;
				}
			}
		}
		
		if(diffbp&&extradata_list!=null){
			for (int i = 0; i < probmatrix2.length; i++) {
				for (int j = 1; j < 4; j++) {
					if (i + j < probmatrix2.length) {
						probmatrix2[i][i + j] = 0;
					}
					if (i - j > 0) {
						probmatrix2[i][i - j] = 0;
					}
				}
			}
		}
		

		System.out.println("Folding...");
		act.setCurrentActivity("Applying grammar");
		
		//starttime = System.currentTimeMillis();
		ResultBundle result = FoldingProject.calcSCFG(act
				.getChildProgress(scfgpart), scfgjobsnr, param.getProb(),
				probmatrix, probmatrix2, executor, verbose, diffbp, entropycalc);
		//System.out.println("Time in SCFG part: " + (System.currentTimeMillis()-starttime));
		result.phyloProbs = probmatrix;
		RNAFoldingTools.writeMatrix(probmatrix, new File("probmatrix.txt"));
		
		//Shut down the executor so we aren't hanging at the end
		executor.shutDown();
		return result;

	}
	
	public static ResultBundle foldMatrix(Progress act, int phylojobsnr,
			int scfgjobsnr, List<char[]> columns,
			List<String> names, double [][] inputMatrix, Parameters param,
			AsynchronousJobExecutor executor, boolean verbose, int execnr,
			List <ExtraData> extradata_list, boolean diffbp, boolean entropycalc)
			throws InterruptedException {
		
		if(columns.size()<2){
			return ResultBundle.tinyBundle();
		}
		
		double scfg_to_phylo_ratio = 0.095 * columns.size()
				/ columns.get(0).length;
		double phylopart = 0.95 / (scfg_to_phylo_ratio + 1);
		if (phylopart < 0 | phylopart > 1) {
			System.err.println("Time estimation resulted "
					+ "in illegal value for phylogenetic calculations: "
					+ phylopart);
			phylopart = 0.475; // in case the estimation results in a
			// "weird number", just set it to 50-50%.
		}
		double scfgpart = 0.95 - phylopart; // max is 95% because processing
		// must happen before and after
		act.setCurrentActivity("Applying evolutionary model");

		
		long starttime = System.currentTimeMillis();
		/*double[][] probmatrix = FoldingProject.createPhyloProb(act
				.getChildProgress(phylopart), phylojobsnr, tree, columns,
				names, columns.size(), param, executor, verbose, execnr);*/

		//If there are extra data, multiply probmatrix with the relevant values
		
		double [][] probmatrix = new double[inputMatrix.length][inputMatrix[0].length];
		for(int i = 0 ; i < probmatrix.length ; i++)
		{
			double sum = 0;
			for(int j = 0 ; j < probmatrix[0].length ; j++)
			{
				probmatrix[i][j] = inputMatrix[i][j];
				sum += probmatrix[i][j];
			}
			//probmatrix[i][i] = 1 - sum;
			System.out.println("XVX"+sum);
		}
		
		double[][] probmatrix2 = null;
		if(diffbp&&extradata_list!=null){
			probmatrix2= new double[probmatrix.length][probmatrix[0].length];
			for (int a = 0; a<probmatrix2.length;a++){
				for(int b = 0;b<probmatrix2[0].length; b++){
					probmatrix2[a][b] = probmatrix[a][b]; 
				}
			}
		}
		
		if(extradata_list!=null){
			System.out.println("Number of auxiliary data items: " + extradata_list.size());
			ArrayList<ForcedConstraints> fcList = new ArrayList<ForcedConstraints>();
			for(ExtraData data:extradata_list){
				if(data instanceof ForcedConstraints){
					//If forcing constraints, collect them into one and remove from the list
					fcList.add((ForcedConstraints) data);
				}
			}
			for(ExtraData data:fcList){
				extradata_list.remove(data);
			}
			if(fcList.size()>0){
				System.out.println("Processing " + fcList.size () + " hard constraint datasets");
				ForcedConstraints fcMerged = new ForcedConstraints();
				fcMerged = fcMerged.combinedForcedConstraints(probmatrix.length, fcList);
				extradata_list.add(fcMerged);
			}

		}
		
		System.out.println("Processing all auxiliary data");
		if(extradata_list!=null){
			for(ExtraData extradata:extradata_list){
				for(int i = 0; i<probmatrix.length; i++){
					for(int j = 0; j<i; j++){
						probmatrix[i][j] *= extradata.getProbabilityGivenInnerPaired(i,j);
						probmatrix[j][i] = probmatrix[i][j];
					}
					probmatrix[i][i] *= extradata.getProbabilityGivenUnpaired(i);
				}
			}
		}
		
		if(diffbp&&extradata_list!=null){
			//If inner and outer basepairs are distinguished
			for(ExtraData extradata:extradata_list){
				for(int i = 0; i<probmatrix2.length; i++){
					for(int j = 0; j<i; j++){
						probmatrix2[i][j] *= extradata.getProbabilityGivenOuterPaired(i,j);
						probmatrix2[j][i] = probmatrix2[i][j];
					}
					probmatrix2[i][i] *= extradata.getProbabilityGivenUnpaired(i);
				}
			}
			
		}
		
		//System.out.println("Time in phylogenetic part: " + (System.currentTimeMillis()-starttime));
		// do not allow basepairs less than 4 nucleotides apart.
		for (int i = 0; i < probmatrix.length; i++) {
			for (int j = 1; j < 4; j++) {
				if (i + j < probmatrix.length) {
					probmatrix[i][i + j] = 0;
				}
				if (i - j > 0) {
					probmatrix[i][i - j] = 0;
				}
			}
		}
		
		if(diffbp&&extradata_list!=null){
			for (int i = 0; i < probmatrix2.length; i++) {
				for (int j = 1; j < 4; j++) {
					if (i + j < probmatrix2.length) {
						probmatrix2[i][i + j] = 0;
					}
					if (i - j > 0) {
						probmatrix2[i][i - j] = 0;
					}
				}
			}
		}
		

		System.out.println("Folding...");
		act.setCurrentActivity("Applying grammar");
		
		//starttime = System.currentTimeMillis();
		ResultBundle result = FoldingProject.calcSCFG(act
				.getChildProgress(scfgpart), scfgjobsnr, param.getProb(),
				probmatrix, probmatrix2, executor, verbose, diffbp, entropycalc);
		//System.out.println("Time in SCFG part: " + (System.currentTimeMillis()-starttime));

		//Shut down the executor so we aren't hanging at the end
		executor.shutDown();
		return result;

	}
	
	public static ResultBundle foldFuzzyAlignment(Progress act, int phylojobsnr,
			int scfgjobsnr, Tree tree, List<FuzzyNucleotide[]> columns, FuzzyAlignment fuzzyAlignment,
			List<String> names, Parameters param,
			AsynchronousJobExecutor executor, boolean verbose, int execnr,
			List <ExtraData> extradata_list, boolean diffbp, boolean entropycalc)
			throws InterruptedException {
		
		if(columns.size()<2){
			return ResultBundle.tinyBundle();
		}
		
		double scfg_to_phylo_ratio = 0.095 * columns.size()
				/ columns.get(0).length;
		double phylopart = 0.95 / (scfg_to_phylo_ratio + 1);
		if (phylopart < 0 | phylopart > 1) {
			System.err.println("Time estimation resulted "
					+ "in illegal value for phylogenetic calculations: "
					+ phylopart);
			phylopart = 0.475; // in case the estimation results in a
			// "weird number", just set it to 50-50%.
		}
		double scfgpart = 0.95 - phylopart; // max is 95% because processing
		// must happen before and after
		act.setCurrentActivity("Applying evolutionary model");

		
		//long starttime = System.currentTimeMillis();
		double[][] probmatrix = FoldingProject.createPhyloProbForFuzzyAlignment(act
				.getChildProgress(phylopart), phylojobsnr, tree, columns, fuzzyAlignment,
				names, columns.size(), param, executor, verbose, execnr);

		//If there are extra data, multiply probmatrix with the relevant values
		
		double[][] probmatrix2 = null;
		if(diffbp&&extradata_list!=null){
			probmatrix2= new double[probmatrix.length][probmatrix[0].length];
			for (int a = 0; a<probmatrix2.length;a++){
				for(int b = 0;b<probmatrix2[0].length; b++){
					probmatrix2[a][b] = probmatrix[a][b]; 
				}
			}
		}
		
		if(extradata_list!=null){
			System.out.println("Number of auxiliary data items: " + extradata_list.size());
			ArrayList<ForcedConstraints> fcList = new ArrayList<ForcedConstraints>();
			for(ExtraData data:extradata_list){
				if(data instanceof ForcedConstraints){
					//If forcing constraints, collect them into one and remove from the list
					fcList.add((ForcedConstraints) data);
				}
			}
			for(ExtraData data:fcList){
				extradata_list.remove(data);
			}
			if(fcList.size()>0){
				System.out.println("Processing " + fcList.size () + " hard constraint datasets");
				ForcedConstraints fcMerged = new ForcedConstraints();
				fcMerged = fcMerged.combinedForcedConstraints(probmatrix.length, fcList);
				extradata_list.add(fcMerged);
			}

		}
		
		System.out.println("Processing all auxiliary data");
		if(extradata_list!=null){
			for(ExtraData extradata:extradata_list){
				for(int i = 0; i<probmatrix.length; i++){
					for(int j = 0; j<i; j++){
						probmatrix[i][j] *= extradata.getProbabilityGivenInnerPaired(i,j);
						probmatrix[j][i] = probmatrix[i][j];
					}
					probmatrix[i][i] *= extradata.getProbabilityGivenUnpaired(i);
				}
			}
		}
		
		if(diffbp&&extradata_list!=null){
			//If inner and outer basepairs are distinguished
			for(ExtraData extradata:extradata_list){
				for(int i = 0; i<probmatrix2.length; i++){
					for(int j = 0; j<i; j++){
						probmatrix2[i][j] *= extradata.getProbabilityGivenOuterPaired(i,j);
						probmatrix2[j][i] = probmatrix2[i][j];
					}
					probmatrix2[i][i] *= extradata.getProbabilityGivenUnpaired(i);
				}
			}
			
		}
		
		//System.out.println("Time in phylogenetic part: " + (System.currentTimeMillis()-starttime));
		// do not allow basepairs less than 4 nucleotides apart.
		for (int i = 0; i < probmatrix.length; i++) {
			for (int j = 1; j < 4; j++) {
				if (i + j < probmatrix.length) {
					probmatrix[i][i + j] = 0;
				}
				if (i - j > 0) {
					probmatrix[i][i - j] = 0;
				}
			}
		}
		
		if(diffbp&&extradata_list!=null){
			for (int i = 0; i < probmatrix2.length; i++) {
				for (int j = 1; j < 4; j++) {
					if (i + j < probmatrix2.length) {
						probmatrix2[i][i + j] = 0;
					}
					if (i - j > 0) {
						probmatrix2[i][i - j] = 0;
					}
				}
			}
		}
		

		System.out.println("Folding...");
		act.setCurrentActivity("Applying grammar");
		
		//starttime = System.currentTimeMillis();
		ResultBundle result = FoldingProject.calcSCFG(act
				.getChildProgress(scfgpart), scfgjobsnr, param.getProb(),
				probmatrix, probmatrix2, executor, verbose, diffbp, entropycalc);
		//System.out.println("Time in SCFG part: " + (System.currentTimeMillis()-starttime));

		//Shut down the executor so we aren't hanging at the end
		executor.shutDown();
		return result;

	}
	
	
	
	private static double[][] createPhyloProbForFuzzyAlignment(Progress act, int userjobsnr,
			Tree tree, List<FuzzyNucleotide[]> columns_char, FuzzyAlignment fuzzyAlignment, List<String> names,
			final int length, Parameters param,
			AsynchronousJobExecutor executor, final boolean verbose, int execnr)
			throws InterruptedException {
		final long starttime = System.nanoTime();
		if(verbose){
		System.out.println("Timer (phylogeny) started. (time: "
				+ (System.nanoTime() - starttime) * 1e-9 + " s)");
		}
		if (verbose) {
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if (verbose) {
			System.out.println("Processing input...");
		}
		if (verbose) {
			act.setCurrentActivity("Evolutionary model: processing input");
		}
		List<FuzzyNucleotide[]> columns = new ArrayList<FuzzyNucleotide[]>();
		
		
		/*
		for (int i = 0; i < columns_char.size(); i++) {
			columns.add(MatrixTools.convertColumn(columns_char.get(i)));
		}*/
		for (int i = 0; i < columns_char.size(); i++) {
			columns.add(columns_char.get(i));
		}
		
		if (verbose) {
			System.out.println("User wish for number of divisions: " + userjobsnr);
		}

		// correct user input
		if (userjobsnr > length) {
			userjobsnr = length;
		} else if (userjobsnr < 1) {
			userjobsnr = 1;
		}
		int nrjobs = userjobsnr;

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if (verbose) {
			System.out
					.println("Actual nr. of divisions in phylogenetic calculations: "
							+ nrjobs);
		}
		if (verbose) {
			System.out.println("Generating matrices for each node...");
		}

		tree.generateLeafList(names); // to speed up calculations later

		final double[][] probmatrix = new double[length][length];
		List<PhyloJobFuzzy> jobs = new ArrayList<PhyloJobFuzzy>();
		// First generate "exp(Rt)" matrix for each node.
		// This is the same for all columns of that node.
		tree.getRoot().calculateChildrenMatrix(param.getSD(), param.getSV(),
				param.getSV1());

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out.println("Generating jobs...");
		}
		if (verbose) {
			act.setCurrentActivity("Evolutionary model: dividing tasks");
		}

		// Create jobs for SINGLE columns
		// now all single columns are in one job
		// Add the only single-column job
		int colcnt = 0;
		PhyloJobFuzzy lastjob = new PhyloJobFuzzy(fuzzyAlignment);
		lastjob.tree = Tree.copyTree(tree); // must copy the entire tree into
		// each phylojob
		lastjob.names = names; // must copy the names into each phylojob (for
		// debugging only)
		lastjob.startcol = 0;
		lastjob.jobid = 0;
		lastjob.endcol = length - 1;
		lastjob.type = false;
		lastjob.param = param;
		for (int col = 0; col < length; col++) {
			// send columns
			lastjob.columns.add(columns.get(col));
			lastjob.columnIndices.add(col);
			colcnt++;
		}
		jobs.add(lastjob);

		// Create jobs for column PAIRS
		// count total number of column pairs
		int paircnt = 0;
		for (int i = 0; i < length; i++) {
			for (int j = i + 1; j < length; j++) {
				paircnt++;
			}
		}
		if(verbose){
			System.out.println("Total number of pairs: " + paircnt);
			System.out.println("Pairs in a job: " + (paircnt) / nrjobs);
		}
		// have to calculate matrices, but only once
		tree.getRoot().calculateChildrenMatrix(param.getDD(), param.getDV(),
				param.getDV1());

		int lastendcol = -1;
		int paircolcnt = 0;
		// add all jobs except last
		for (int jobnr = 0; jobnr < nrjobs - 1; jobnr++) {
			PhyloJobFuzzy job = new PhyloJobFuzzy(fuzzyAlignment);
			job.tree = Tree.copyTree(tree); // must copy the entire tree into
			// each phylojob
			job.names = names; // must copy the names into each phylojob
			job.type = true;
			job.param = param;
			job.jobid = jobnr + 1;

			int col1 = lastendcol + 1;
			job.startcol = col1;
			paircolcnt += length - col1;
			col1--;
			while (paircolcnt < ((long)(jobnr + 1) * (paircnt)) / nrjobs) {
				col1++;
				if (col1 == length) { // no more columns left, finish
					col1--;
					break;
				}
				job.columns.add(columns.get(col1));
				job.columnIndices.add(col1);
				paircolcnt += length - col1;
			}

			job.endcol = col1;
			lastendcol = col1;
			for (int col2 = job.startcol + 1; col2 < length; col2++) {
				// send pairing columns
				job.columns2.add(columns.get(col2));
				job.columnIndices2.add(col2);
			}
			jobs.add(job);
			paircolcnt -= length - col1;
		}
		// Add the last job
		lastjob = new PhyloJobFuzzy(fuzzyAlignment);
		lastjob.tree = Tree.copyTree(tree); // must copy the entire tree into
		// each phylojob
		lastjob.names = names; // must copy the names into each phylojob
		lastjob.startcol = lastendcol + 1;
		lastjob.endcol = length;
		lastjob.type = true;
		lastjob.param = param;
		lastjob.jobid = nrjobs;
		for (int col1 = lastjob.startcol; col1 < length; col1++) {
			// send columns
			lastjob.columns.add(columns.get(col1));
			lastjob.columnIndices.add(col1);
		}
		for (int col2 = lastjob.startcol + 1; col2 < length; col2++) {
			// send pairing columns
			lastjob.columns2.add(columns.get(col2));
			lastjob.columnIndices2.add(col2);
		}
		jobs.add(lastjob);

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out
					.println("Total number of jobs in phylogenetic calculations: "
							+ jobs.size());
		}
		if (verbose) {
			System.out.println("Executing jobs...");
		}

		// start executing the jobs
		// counts how many phylojobs are done		
		final AtomicInteger finishedphylojobscount = new AtomicInteger(0); 

		act.setCurrentActivity("Evolutionary model: " +
					"calculating column probabilities");
		
		long gridstarttime = System.nanoTime();
		// execute single column jobs
		Progress singleColAct = act.getChildProgress(0.1);
		for (int jobnr = 0; jobnr < 1; jobnr++) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			PhyloJobFuzzy job = jobs.get(jobnr);
			//final int jn = jobnr;
			final int startcol = job.startcol;

			final Progress jobAct = singleColAct
					.getChildProgress((float) job.columns.size()
							/ (float) colcnt);
			executor.startExecution(job, new JobListener() {
				public void jobFinished(JobResults result) {
				} // doesn't happen here

				public void jobFinished(List<ResultBundle> result) {
				} // doesn't happen here

				public void jobFinished(double[][] result) {
					for (int col = 0; col < result.length; col++) {
						probmatrix[col + startcol][col + startcol] = result[col][0];
					}
					finishedphylojobscount.incrementAndGet();
					jobAct.setProgress(1.0);
				}
			});
		}

		Progress doubleColAct = act.getChildProgress(0.9);

		// execute double column jobs
		for (int jobnr = 1; jobnr < jobs.size(); jobnr++) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			PhyloJobFuzzy job = jobs.get(jobnr);
			//final int jn = jobnr;
			final int startcol = job.startcol;
			final int endcol = job.endcol;
			final Progress jobAct = doubleColAct
					.getChildProgress((float) job.columns.size()
							/ (float) length);
			if (job.columns2.size() != 0) {
				executor.startExecution(job, new JobListener() {
					public void jobFinished(JobResults result) {
					} // doesn't happen here

					public void jobFinished(List<ResultBundle> result) {
					} // doesn't happen here

					public void jobFinished(double[][] result) {
						for (int col1 = startcol; col1 < endcol + 1; col1++) {
							for (int col2 = col1 + 1; col2 < length; col2++) {
								probmatrix[col1][col2] = result[col1 - startcol][col2
										- startcol - 1];
								probmatrix[col2][col1] = probmatrix[col1][col2]; // make
								// it
								// symmetric
							}
						}
						finishedphylojobscount.incrementAndGet();
						jobAct.setProgress(1.0);
					}
				});
			} else {
				finishedphylojobscount.incrementAndGet();
				jobAct.setProgress(1.0);
			}
		}

		// wait for last thread to finish
		// wait for jobs to finish
		while (finishedphylojobscount.get() < jobs.size()) {
			Thread.sleep(100);
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
		}
		singleColAct.setProgress(1.0);
		doubleColAct.setProgress(1.0);
		act.setProgress(1.0);
		//act.setCurrentActivity("Evolutionary model: done");

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		System.out.println("TOTAL TIME ELAPSED IN PHYLOGENETIC PART: "
				+ (int)((System.nanoTime() - starttime) * 1e-9) + " seconds ");
		
		if(verbose){
			System.out.println("                    ...of which distributed: "
				+ (System.nanoTime() - gridstarttime) * 1e-9 + " seconds");
		}

		//System.out.println();

		// Result contains the a priori probability distribution matrix.
		return probmatrix;
	}

	private static ResultBundle calcSCFG(Progress act, int userjobsnr,
			double[][] prob, double[][] probmatrix, double[][] probmatrix2,
			AsynchronousJobExecutor executor, final boolean verbose, final boolean diffbp, final boolean entropycalc)
			throws InterruptedException {
		PointRes tmp = new PointRes(0, 0);
		final long starttime = System.nanoTime();
		

		
		double entropyVal = 0;
		double entropyPercOfMax = 0;
		double entropyMax = 0;

		if (verbose) {
			System.out.println("Timer (SCFG) started. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
			System.out.println("Processing input...");
		}

		act.setCurrentActivity("Applying grammar: processing input");
		

		int length = probmatrix.length;

		// First corrections of user input
		if (userjobsnr > length) {
			userjobsnr = length - 1;
		} else if (userjobsnr < 1) {
			userjobsnr = 1;
		} else {
		}

		// Now calculate number of jobs
		int distance = length / userjobsnr + 1; // the +1 is to prevent
		// generation of superfluous
		// jobs

		int firstrowspaces = length + 1;

		int isthereextra = (firstrowspaces % distance != 0) ? 1 : 0; 
		// triangle is compelete iff isthereextra = 0.
		int nrdivisions = ((firstrowspaces - firstrowspaces % distance) / distance)
				+ isthereextra;
		// nrdivisions = how many jobs should there be in 1st row?
		int nrsectors = nrdivisions * (nrdivisions + 1) / 2; 
		// total number of jobs in full triangle

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out.println("User wish for number of divisions: "
					+ userjobsnr);
			System.out
					.println("Number of sectors in inside-outside calculations: "
							+ nrsectors);
			System.out
					.println("Distance for inside-outside has been calculated as: "
							+ distance);
			System.out
					.println("Corrected number of jobs in first row of inside-outside: "
							+ nrdivisions);
			System.out
					.println("Total number of inside-outside jobs to generate: "
							+ nrsectors);
		}

		int interestingpoints = length * length / 2;
		int totalpoints = nrsectors * distance * distance;
		int extrapoints = totalpoints - interestingpoints;
		float fraction = (float) extrapoints * 100 / interestingpoints;
		if(verbose){
			System.out.println("Divisions = " + nrdivisions
				+ ", Interesting points: " + interestingpoints
				+ ", Extra points: " + extrapoints + ", Total points: "
				+ totalpoints + ", " + "Fractional extra: " + fraction + "%");
		}
		if (verbose) {
			System.out.println("Generating sectors: ");
		}
		act.setCurrentActivity("Applying grammar: generating sectors");

		Sector top = SectorGenerator.GenerateSectors(nrsectors, distance,
				nrdivisions, length,diffbp);

		final Master master = new Master(top, prob);

		// sector: set basepairs (phylogenetic probabilities) for inside sectors
		PointRes number = new PointRes(0, 0);
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < length - i; j++) {
				Sector sectoro = findSector(i, j, master.bottom);
				int so = findPointST(i, j, sectoro, distance)[0];
				int to = findPointST(i, j, sectoro, distance)[1];
				number.setToDouble(probmatrix[i][i + j]);
				sectoro.setBasePairs(so, to, number);
			}
		}
		
		if(diffbp&&probmatrix2!=null){
			//if differentiating inner and outerbasepairs, create a new bp matrix per sector
			for (int i = 0; i < length; i++) {
				for (int j = 0; j < length - i; j++) {
					Sector sectoro = findSector(i, j, master.bottom);
					int so = findPointST(i, j, sectoro, distance)[0];
					int to = findPointST(i, j, sectoro, distance)[1];
					number.setToDouble(probmatrix2[i][i + j]);
					sectoro.setBasePairs2(so, to, number);
				}
			}
		}
		
		
		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}


		act.setCurrentActivity("Applying grammar: inside algorithm");
		
		if (verbose) {
			System.out.println("Calculating inside values... ");
		}

		final AtomicInteger finishedinsidejobscount = new AtomicInteger(0); // counts
		// how many inside jobs are done
		master.CreateInsideJobChannel();
		Progress insideAct = act.getChildProgress(0.32);
		//final long gridstarttime = System.nanoTime();
		while (master.unProcessedInsideSectors()) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			CYKJob cYKJob = master.takeNextInsideJob(); // This call will block
			// until a job is ready
			final Progress jobAct = insideAct.getChildProgress(1.0 / nrsectors);
			// System.out.println(cYKJob.sectorid + " " + " 0 " +
			// (System.nanoTime()-starttime));
			final int sectorNumber = cYKJob.getSectorid();
			executor.startExecution(cYKJob, new JobListener() {
				public void jobFinished(JobResults result) {
					master.setInsideResult(sectorNumber, result);
					jobAct.setProgress(1.0);
					// System.out.println(sectorNumber + " " + " 1 " +
					// (System.nanoTime()-starttime));
					finishedinsidejobscount.incrementAndGet();
				}

				public void jobFinished(double[][] result) {
				}// doesn't happen here

				public void jobFinished(List<ResultBundle> result) {
				} // doesn't happen here
			});
		}

		// wait for last job to finish
		while (finishedinsidejobscount.get() < nrsectors) {
			Thread.sleep(100);
			if(act.shouldStop()){	
				executor.shutDown();
			}
			act.checkStop();
		}

		insideAct.setProgress(1.0);

		if(verbose){
			System.out.println("Top inside: "
				+ master.top.getInsideMatrixS().getProb(distance - 1,
						length - 1 - master.top.pos[1]));
		}
		//long insidetime = System.nanoTime();
		//System.out.print(nrdivisions + " - ");
		//System.out.println("TOTAL TIME ELAPSED IN INSIDE PART (ALL): "
		//		+ (insidetime - starttime) * 1e-9 + " seconds ");
		//double gridtime1 = (insidetime - gridstarttime) * 1e-9;
		//System.out.println("TOTAL TIME ELAPSED IN INSIDE PART (DISTRIBUTED): "
		//		+ gridtime1 + " seconds ");

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		act.setCurrentActivity("Applying grammar: outside algorithm");
		
		if (verbose) {
			System.out.println("Calculating outside values...");
		}
		master.CreateOutsideJobChannel();
		final AtomicInteger finishedoutsidejobscount = new AtomicInteger(0); // counts
		// how many outside jobs are done

		// outside algorithm
		// set basepairs for outside algorithm (they are shifted relative to
		// inside algo)
		number.setToFloat(0);
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < length - i; j++) {
				Sector sectoro = findSector(i, j, master.bottom);
				int so = findPointST(i, j, sectoro, distance)[0];
				int to = findPointST(i, j, sectoro, distance)[1];
				if (i > 0 && (j + i + 1) < length) {
					number.setToDouble(probmatrix[i - 1][i + j + 1]);
				} else if (i == 0) {
					number.setToFloat(0);
				} else if (j + i + 1 == length) {
					number.setToFloat(0);
				}
				sectoro.setBasePairs(so, to, number);
			}
		}
		
		if(diffbp){
			for (int i = 0; i < length; i++) {
				for (int j = 0; j < length - i; j++) {
					Sector sectoro = findSector(i, j, master.bottom);
					int so = findPointST(i, j, sectoro, distance)[0];
					int to = findPointST(i, j, sectoro, distance)[1];
					if (i > 0 && (j + i + 1) < length) {
						number.setToDouble(probmatrix2[i - 1][i + j + 1]);
					} else if (i == 0) {
						number.setToFloat(0);
					} else if (j + i + 1 == length) {
						number.setToFloat(0);
					}
					sectoro.setBasePairs2(so, to, number);
				}
			}
		}
		
		//final long outsidegridstarttime = System.nanoTime();
		Progress outsideAct = act.getChildProgress(0.50);
		while (master.unProcessedOutsideSectors()) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			final Progress jobAct = outsideAct
					.getChildProgress(1.0 / nrsectors);
			CYKJob cYKJob = master.takeNextOutsideJob(); // This call will block
			// until a job is
			// ready
			final int sectorNumber = cYKJob.getSectorid();
			executor.startExecution(cYKJob, new JobListener() {
				// System.out.println(sectorNumber + " " + " 2 " +
				// (System.nanoTime()-starttime));
				public void jobFinished(List<ResultBundle> result) {
				} // doesn't happen here

				public void jobFinished(JobResults result) {
					master.setOutsideResult(sectorNumber, result);
					jobAct.setProgress(1.0);
					// System.out.println(sectorNumber + " " + " 3 " +
					// (System.nanoTime()-starttime));
					finishedoutsidejobscount.incrementAndGet();
				}

				public void jobFinished(double[][] result) {
				}// doesn't happen here
			});
		}

		// wait for last job to finish
		while (finishedoutsidejobscount.get() < nrsectors) {
			Thread.sleep(100);
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
		}
		outsideAct.setProgress(1.0);

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

	//	long outsidetime = (System.nanoTime());
	//	System.out.print(nrdivisions + " - ");
	//	System.out.println("TOTAL TIME ELAPSED IN OUTSIDE PART: "
	//			+ (outsidetime - insidetime) * 1e-9 + " seconds ");
	//	double gridtime2 = (outsidetime - outsidegridstarttime) * 1e-9;
	//	System.out.println("TOTAL TIME ELAPSED IN OUTSIDE PART (DISTRIBUTED): "
	//			+ gridtime2 + " seconds ");

		if (verbose) {
			System.out.println("Setting basepairs...");
		}

		act.setCurrentActivity("Applying grammar: setting basepairs");

		number.setToFloat(0);
		PointRes number2 = new PointRes(0, 0);

		//Expected rule frequencies
		PointRes EfLdFd = new PointRes(0,0);
		PointRes EfFdFd = new PointRes(0,0);
		PointRes EfSL = new PointRes(0,0);
		PointRes EfSLS = new PointRes(0,0);
		PointRes EfLs = new PointRes(0,0);
		PointRes EfFLS = new PointRes(0,0);
		//Expected nonterminal frequencies (for probability reestimation)
		PointRes EfS = new PointRes(0,0);
		PointRes EfL = new PointRes(0,0);
		PointRes EfF = new PointRes(0,0);
		
		PointRes topInsideS = new PointRes(0,0);
		topInsideS.copyFrom(master.top.getInsideMatrixS().fetchProb(
				distance - 1, length - 1 - master.top.pos[1], tmp));
		
		PointRes insideval = new PointRes(0,0);
		PointRes outsideval = new PointRes(0,0);
		
		PointRes entropyLikelihood = new PointRes(0,0); 
		
		//MatrixTools.print(probmatrix);
		// calculate basepair probabilities & basepair rule entropies
		for (int i = 0; i < length; i++) {
			for (int j = 1; j < length - i; j++) {
				Sector sectoro = findSector(i, j, master.bottom);
				Sector sectori = findSector(i + 1, j - 2, master.bottom);
				int so = findPointST(i, j, sectoro, distance)[0];
				int to = findPointST(i, j, sectoro, distance)[1];
				int si = findPointST(i + 1, j - 2, sectori, distance)[0];
				int ti = findPointST(i + 1, j - 2, sectori, distance)[1];
				
				//L->dFd
				number.copyFrom(sectoro.getOutsideMatrixL().fetchProb(so, to,
						tmp));
				number.multiply(sectori.getInsideMatrixF().fetchProb(si, ti,
						tmp), prob[1][0]);
				if(diffbp){
					number.multiply(probmatrix2[i][i+j]);
				}else{
					number.multiply(probmatrix[i][i+j]);
				}
				EfLdFd.add(number);
				number.divide(topInsideS);
				
				//F->dFd
				number2.copyFrom(sectoro.getOutsideMatrixF().fetchProb(so, to,
						tmp));
				number2.multiply(sectori.getInsideMatrixF().fetchProb(si, ti,
						tmp), prob[2][0]);
				number2.multiply(probmatrix[i][i+j]);
				EfFdFd.add(number2);
				number2.divide(topInsideS);
				
				//Sum of expectation of the two rules is expectation of basepairing rule
				number.add(number2);
				
				sectoro.setBasePairs(so, to, number);
				if(number.toFloat()>1){
					if(verbose){
						System.err.println("Warning: Pd(" + i + ", " + j +
							") = " + number + " > 1! (Using 1...)" );
					}
					number.setToFloat(1);
				}
				if(probmatrix[i][i+j]!=0){
					number.multiply(log2(1/probmatrix[i][i+j]));
					//System.out.println(i + ", " + (i+j) + " pair entropy: " + number);
					entropyLikelihood.add(number);
				}
				
				//System.out.print(number.toFloat() + " ");
				//System.out.format("%16.8e", number.toDouble());
			}
		}
		
		
		
	
		if (verbose) { 
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}


		if (verbose) {
			System.out.println("Calculating single base probabilities...");
		}

		// calculate single basepairing
		final PointRes one = new PointRes(1);
		final PointRes minusone = new PointRes(-1);
		final PointRes zero = new PointRes(0);

		for (int i = 0; i < length; i++) {
			PointRes singlebaseprobpoint = new PointRes(0);
			for (int k = 1; k < length - i; k++) { // pairs to the right
				Sector basepairsector = findSector(i, k, master.bottom);
				int[] basepairpoint = findPointST(i, k, basepairsector,
						distance);
				singlebaseprobpoint.add(basepairsector.getMatrixBpVal(
						basepairpoint[0], basepairpoint[1], tmp));
			}

			for (int k = 0; k < i; k++) { // pairs to the left
				Sector basepairsector = findSector(k, i - k, master.bottom);
				int[] basepairpoint = findPointST(k, i - k, basepairsector,
						distance);
				singlebaseprobpoint.add(basepairsector.getMatrixBpVal(
						basepairpoint[0], basepairpoint[1], tmp));
			}

			tmp.copyFrom(minusone);
			singlebaseprobpoint.multiply(tmp);
			tmp.copyFrom(one);
			tmp.add(singlebaseprobpoint);
			singlebaseprobpoint.copyFrom(tmp);

			if (singlebaseprobpoint.isSignificantlyLessThanZero()) {
				System.err.println("Illegal Ps(" + i + ") = "
						+ singlebaseprobpoint
						+ ": (significantly less than zero); using 0");
				singlebaseprobpoint.setToFloat(0);
			}

			Sector sec = findSector(i, 0, master.bottom);
			int[] sbasepoint = findPointST(i, 0, sec, distance);
			sec.setBasePairs(sbasepoint[0], sbasepoint[1], singlebaseprobpoint);
			if(probmatrix[i][i]!=0){
				singlebaseprobpoint.multiply(log2(1/probmatrix[i][i]));
				entropyLikelihood.add(singlebaseprobpoint);
			}
			//System.out.format("%16.8e", singlebaseprobpoint.toDouble());
		}

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		
		if(entropycalc){
			for (int i = 0; i < length; i++) {
				//Entropy for L->s
				Sector sector = findSector(i,0, master.bottom);
				int s = findPointST(i,0,sector,distance)[0];
				int t = findPointST(i,0,sector,distance)[1];
				number.setToDouble(probmatrix[i][i]);
				outsideval.copyFrom(sector.getOutsideMatrixL().fetchProb(s, t, tmp));
				number.multiply(outsideval);
				number.multiply(prob[1][1]);
				//number.divide(topInsideS);
				EfLs.add(number);
				
				for (int j = 0; j < length - i; j++) {
					Sector sectoro = findSector(i, j, master.bottom);
					int so = findPointST(i, j, sectoro, distance)[0];
					int to = findPointST(i, j, sectoro, distance)[1];
					
					insideval.copyFrom(sectoro.getInsideMatrixS().fetchProb(so,to,tmp));
					outsideval.copyFrom(sectoro.getOutsideMatrixS().fetchProb(so, to,tmp));
					outsideval.multiply(insideval);
					EfS.add(outsideval);
					insideval.copyFrom(sectoro.getInsideMatrixL().fetchProb(so,to,tmp));
					outsideval.copyFrom(sectoro.getOutsideMatrixL().fetchProb(so, to,tmp));
					outsideval.multiply(insideval);
					EfL.add(outsideval);
					insideval.copyFrom(sectoro.getInsideMatrixF().fetchProb(so,to,tmp));
					outsideval.copyFrom(sectoro.getOutsideMatrixF().fetchProb(so, to,tmp));
					outsideval.multiply(insideval);
					EfF.add(outsideval);
					
					insideval.copyFrom(sectoro.getInsideMatrixL().fetchProb(so,to,tmp));
					outsideval.copyFrom(sectoro.getOutsideMatrixS().fetchProb(so, to,tmp));
					outsideval.multiply(insideval);
					outsideval.multiply(prob[0][1]);
					//outsideval.divide(topInsideS);
					EfSL.add(outsideval);
				}
			}

// calculate bifurcation probabilities NOTE O(n^3) complexity		
	/*PointRes insidevalue1 = new PointRes(0,0);
		PointRes insidevalue2 = new PointRes(0,0);
		for (int i = 0; i < length; i++) {	
			for (int j = 0; j < length - i; j++) {
				for(int k = 0; k<j; k++){
					//i+k is between i and i+j.
					Sector sectoro = findSector(i, j, master.bottom);
					Sector sectori1 = findSector(i, k, master.bottom);
					Sector sectori2 = findSector(i+k+1, j-k-1, master.bottom);
					//System.out.println(i + " " + (i+k) + " " + (i+j));
					//Outside region (i,i+j)  
					int so = findPointST(i, j, sectoro, distance)[0];
					int to = findPointST(i, j, sectoro, distance)[1];
					//Inside region (i,i+k)
					int si1 = findPointST(i, k, sectori1, distance)[0];
					int ti1 = findPointST(i, k, sectori1, distance)[1];
					//Inside region (i+k+1, j)
					int si2 = findPointST(i+k+1, j-k-1, sectori2, distance)[0];
					int ti2 = findPointST(i+k+1, j-k-1, sectori2, distance)[1];
					//Common to both bifurcation rules 
					insidevalue1.copyFrom(sectori1.getInsideMatrixL().fetchProb(si1,ti1,tmp));
					insidevalue2.copyFrom(sectori2.getInsideMatrixS().fetchProb(si2,ti2,tmp));					
					//S->LS: prob[0][0]
					number.copyFrom(sectoro.getOutsideMatrixS().fetchProb(so, to,	tmp));
					number.multiply(insidevalue1);
					number.multiply(insidevalue2);
					number.multiply(prob[0][0]);
					//number.divide(topInsideS);
					EfSLS.add(number);
					//F->LS: prob[2][1]
					number.copyFrom(sectoro.getOutsideMatrixF().fetchProb(so, to,	tmp));
					number.multiply(insidevalue1);
					number.multiply(insidevalue2);
					number.multiply(prob[2][1]);
					//number.divide(topInsideS);
					EfFLS.add(number);
				}
			}
		}
*/
		//This is to save O(n^3) execution time
		EfSLS.copyFrom(EfS);
		EfSLS.subtract(EfSL,tmp);
		EfFLS.copyFrom(EfF);
		EfFLS.subtract(EfFdFd,tmp);
		
		EfS.divide(topInsideS);
		EfL.divide(topInsideS);
		EfF.divide(topInsideS);
		EfSLS.divide(topInsideS);
		EfSL.divide(topInsideS);
		EfLdFd.divide(topInsideS);
		EfLs.divide(topInsideS);
		EfFLS.divide(topInsideS);
		EfFdFd.divide(topInsideS);
		
		PointRes entropySLS = new PointRes(EfSLS);
		PointRes entropySL = new PointRes(EfSL);
		PointRes entropyLdFd = new PointRes(EfLdFd);
		PointRes entropyLs = new PointRes(EfLs);
		PointRes entropyFdFd = new PointRes(EfFdFd);
		PointRes entropyFLS = new PointRes(EfFLS);
			
		if(verbose){
		System.out.println("Expected nonterminal freq S = " + EfS.toDouble() );
		System.out.println("Expected nonterminal freq L = " + EfL.toDouble() );
		System.out.println("Expected nonterminal freq F = " + EfF.toDouble() );
			
		System.out.println("Expected rule freq S->LS = " + EfSLS.toDouble() );
		System.out.println("Expected rule freq S->L = " + EfSL.toDouble() );
		System.out.println("Expected rule freq L->dFd = " + EfLdFd.toDouble() );
		System.out.println("Expected rule freq L->s = " + EfLs.toDouble() );
		System.out.println("Expected rule freq F->LS = " + EfFLS.toDouble() );
		System.out.println("Expected rule freq F->dFd = " + EfFdFd.toDouble() );
		}

		//EfSLS.divide(EfS);
		EfSL.divide(EfS);
		EfSLS.setToDouble(1-EfSL.toDouble());
		EfLdFd.divide(EfL);
		EfLs.divide(EfL);
		EfFdFd.divide(EfF);
		//EfFLS.divide(EfF);
		EfFLS.setToDouble(1-EfFdFd.toDouble());
		
		if(verbose){
		System.out.println("Reestimation prob S->LS = " + EfSLS.toDouble() );
		System.out.println("Reestimation prob S->L = " + EfSL.toDouble() );
		System.out.println("Reestimation prob L->dFd = " + EfLdFd.toDouble() );
		System.out.println("Reestimation prob L->s = " + EfLs.toDouble() );
		System.out.println("Reestimation prob F->LS = " + EfFLS.toDouble() );
		System.out.println("Reestimation prob F->dFd = " + EfFdFd.toDouble() );
		}
		/*
		System.out.println("CHECK: Sum of rules from nonterminals:");
		System.out.println("S: " + (float)( EfSLS.toDouble() + EfSL.toDouble() ));
		System.out.println("L: " + (float)( EfLdFd.toDouble() + EfLs.toDouble() ));
		System.out.println("F: " + (float)( EfFdFd.toDouble() + EfFLS.toDouble() ));
		*/
				
		entropySLS.multiply(new PointRes(log2(1/prob[0][0])));//.divide(LogTwo);
		entropySL.multiply(new PointRes(log2(1/prob[0][1])));//.divide(LogTwo);
		entropyLs.multiply(new PointRes(log2(1/prob[1][1])));//.divide(LogTwo);
		entropyFLS.multiply(new PointRes(log2(1/prob[2][1])));//.divide(LogTwo);
		entropyLdFd.multiply(new PointRes(log2(1/prob[1][0])));//.divide(LogTwo);
		entropyFdFd.multiply(new PointRes(log2(1/prob[2][0])));//.divide(LogTwo);
	
		PointRes entropy = new PointRes(entropySLS);
		entropy.add(entropySL);
		entropy.add(entropyLdFd);
		entropy.add(entropyLs);
		entropy.add(entropyFLS);
		entropy.add(entropyFdFd);
		if(verbose){System.out.println("Rule entropy component: " + entropy);}

		PointRes logTopInsideS = new PointRes(topInsideS);
		logTopInsideS.takeLog2();
		//logTopInsideS.multiply(-1);
		entropy.add(logTopInsideS);
		if(verbose){System.out.println("Log top inside S component: " + logTopInsideS);}
		//System.out.println("Log top inside S + rule entropy: " + entropy);
		if(verbose){System.out.println("Likelihood entropy component: " + entropyLikelihood);}
		entropy.add(entropyLikelihood);
		
		System.out.println("ENTROPY: " + entropy + " = " + entropy.toDouble());
		
		
		DecimalFormat df2 = new DecimalFormat( "#########0.00" );
		double maxentropy_nr = 0.142- 1.5*Math.log(length)/Math.log(2) + 1.388*length;
		double maxentropy = new Double(df2.format(maxentropy_nr)).doubleValue();
		double percent_nr = (entropy.toDouble()/maxentropy_nr)*100;
		double percent = new Double(df2.format(percent_nr)).doubleValue();
		
		entropyVal = entropy.toDouble();
		entropyPercOfMax = percent;
		entropyMax = maxentropy;
		
		System.out.println("(which is " + percent + "% of the maximum entropy, " + maxentropy + ")");
		
		}
		
		if (verbose) {
			System.out.println("Cleaning memory...");
		}

		// assign null's to all inside-outside to enforce clearing of memory in
		// case it didn't happen
		Sector thissec = master.bottom;
		thissec.clearAllInside();
		thissec.clearAllOutside();
		while (thissec.next != null) {
			thissec = thissec.next;
			thissec.clearAllInside();
			thissec.clearAllOutside();
		}

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}		
		
		if (verbose) {
			System.out.println("Calculating expectation values...");
		}

			act.setCurrentActivity("Applying grammar: setting expectation values");

		// calculate expectation values
		master.CreateExpectationJobChannel();
		final AtomicInteger finishedexpectationjobscount = new AtomicInteger(0); // counts
		// how many exp jobs are done
		Progress expectAct = act.getChildProgress(0.11);
		//final long expgridstarttime = System.nanoTime();
		while (master.unProcessedExpectationSectors()) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			final Progress jobAct = expectAct.getChildProgress(1.0 / nrsectors);
			CYKJob cYKJob = master.takeNextExpectationJob(); // This call will
			// block until a job is ready 
			// System.out.println(cYKJob.sectorid + " " + " 4 " +
			// (System.nanoTime()-starttime));
			final int sectorNumber = cYKJob.getSectorid();
			executor.startExecution(cYKJob, new JobListener() {
				public void jobFinished(JobResults result) {
					master.setExpectationResult(sectorNumber, result);
					jobAct.setProgress(1.0);
					// System.out.println(sectorNumber + " " + " 5 " +
					// (System.nanoTime()-starttime));
					finishedexpectationjobscount.incrementAndGet();
				}

				public void jobFinished(double[][] result) {
				}// doesn't happen here

				public void jobFinished(List<ResultBundle> result) {
				} // doesn't happen here
			});
		}
		// wait for last job to finish
		while (finishedexpectationjobscount.get() < nrsectors) {
			Thread.sleep(100);
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
		}
		expectAct.setProgress(1.0);

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if(verbose){
		System.out.println("Top expectation: "
				+ master.top.getExpectationMatrix().getProb(distance - 1,
						length - 1 - master.top.pos[1]));
		}
		// System.out.print(nrdivisions + " - ");
		//System.out.println("TOTAL TIME ELAPSED IN EXPECTATION PART: "
		//		+ (System.nanoTime() - outsidetime) * 1e-9 + " seconds ");
		//double gridtime3 = (System.nanoTime() - expgridstarttime) * 1e-9;
		//System.out
		//		.println("TOTAL TIME ELAPSED IN EXPECTATION PART (DISTRIBUTED): "
		//				+ gridtime3 + " seconds ");

		if (verbose) {
			System.out.println("Backtracking to find structure...");
		}

		act.setCurrentActivity("Applying grammar: backtracking");
		

		char[] structure = new char[length];
		float[] reliability = new float[length];
		LinkedList<int[]> checkPoint = new LinkedList<int[]>(); 

		int[] startpoint = new int[2];
		startpoint[0] = 0;
		startpoint[1] = length - 1;
		checkPoint.addFirst(startpoint);

		PointRes expectationvalue = new PointRes(0, 0);
		PointRes otherexpectationvalue = new PointRes(0, 0);
		PointRes thirdexpectationvalue = new PointRes(0, 0);
		PointRes basepairvalue = new PointRes(0, 0);
		PointRes tmp2 = new PointRes(0, 0);
		
		//int[] pairing = new int[length]; //this creates the pairing array
		//for(int i = 0; i<length; i++){
		//	pairing[i] = -1;
		//}
		
		while (!checkPoint.isEmpty()) {
			// Get the first element from the LinkedList.
			int[] point = checkPoint.removeFirst();
			// Check if this point is a single base.
			Sector sec = findSector(point[0], point[1] - point[0],
					master.bottom);
			int[] stpoint = findPointST(point[0], point[1] - point[0], sec,
					distance);
			expectationvalue.copyFrom(sec.getMatrixExpectationTpVal(stpoint[0],
					stpoint[1], tmp));
			basepairvalue.copyFrom(sec.getMatrixBpVal(stpoint[0], stpoint[1],
					tmp));
			if (point[0] == point[1]) {
				structure[point[0]] = '.';
				reliability[point[0]] = basepairvalue.toFloat();
				continue;
			}
			// Check if this point is paired.
			Sector othersec = findSector(point[0] + 1, point[1] - 1
					- (point[0] + 1), master.bottom);
			int[] otherstpoint = findPointST(point[0] + 1, point[1] - 1
					- (point[0] + 1), othersec, distance);
			otherexpectationvalue.copyFrom(othersec.getMatrixExpectationTpVal(
					otherstpoint[0], otherstpoint[1], tmp));

			tmp.copyFrom(basepairvalue);
			tmp.multiply(2);
			tmp2.copyFrom(otherexpectationvalue);
			tmp2.add(tmp);

			if (expectationvalue.equals(tmp2)) {
				// do the pairing
				structure[point[0]] = '(';
				structure[point[1]] = ')';
				//pairing[point[0]] = point[1];
				//pairing[point[1]] = point[0]; 
				
				reliability[point[0]] = basepairvalue.toFloat();
				reliability[point[1]] = basepairvalue.toFloat();

				// push next point in the list
				int[] nextpoint = new int[2];
				nextpoint[0] = point[0] + 1;
				nextpoint[1] = point[1] - 1;
				checkPoint.addFirst(nextpoint);
				continue;
			} else { // if it isn't basepaired then it is bifuricated; find the
				// two points that gave the value.
				for (int k = point[0]; k < point[1]; k++) { // "k from i to j-1"
					othersec = findSector(point[0], k - point[0], master.bottom);
					otherstpoint = findPointST(point[0], k - point[0],
							othersec, distance);
					otherexpectationvalue.copyFrom(othersec
							.getMatrixExpectationTpVal(otherstpoint[0],
									otherstpoint[1], tmp));

					Sector thirdsec = findSector(k + 1, point[1] - (k + 1),
							master.bottom);
					int[] thirdpoint = findPointST(k + 1, point[1] - (k + 1),
							thirdsec, distance);
					thirdexpectationvalue.copyFrom(thirdsec
							.getMatrixExpectationTpVal(thirdpoint[0],
									thirdpoint[1], tmp));
					tmp.copyFrom(otherexpectationvalue);
					tmp.add(thirdexpectationvalue);
					if (expectationvalue.equals(tmp)) {
						int[] nextpoint1 = new int[2];
						nextpoint1[0] = point[0];
						nextpoint1[1] = k;

						int[] nextpoint2 = new int[2];
						nextpoint2[0] = k + 1;
						nextpoint2[1] = point[1];

						checkPoint.addFirst(nextpoint2);
						checkPoint.addFirst(nextpoint1);
						break;
					}
				}
				continue;
			}
			
		}
		
		//System.out.println();
		//for(int i = 0; i<length; i++){
		//	System.out.format("%08d", pairing[i]);
		//}

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out.println("Finalizing results... ");
		}

		act.setCurrentActivity("Applying grammar: finalizing");
		

		// copy basepairs to square matrix
		float[][] expectationvalues = new float[length][length];
		float[][] basepairprob = new float[length][length];
		float[] singlebaseprob = new float[length];

		for (int j = 1; j < length; j++) {
			for (int i = 0; i < length - j; i++) {
				Sector sectoro = findSector(i, j, master.bottom);
				int so = sectoro.pos[0] + distance - i - 1;
				int to = j + i - (sectoro.pos[1] + sectoro.pos[0]);
				tmp.copyFrom(sectoro.getMatrixBpVal(so, to, tmp2));
				basepairprob[i][i + j] = tmp.toFloat();
				basepairprob[i + j][i] = basepairprob[i][i + j];
				// make it symmetric
				tmp.copyFrom(sectoro.getMatrixExpectationTpVal(so, to, tmp2));
				expectationvalues[i][i + j] = tmp.toFloat();
				expectationvalues[i + j][i] = expectationvalues[i][i + j]; 
				// make it symmetric
			}
		}
		for (int i = 0; i < length; i++) {
			Sector sec = findSector(i, 0, master.bottom);
			int s = sec.pos[0] + distance - i - 1;
			int t = i - (sec.pos[1] + sec.pos[0]);
			tmp.copyFrom(sec.getMatrixBpVal(s, t, tmp2));
			singlebaseprob[i] = tmp.toFloat();
			tmp.copyFrom(sec.getMatrixExpectationTpVal(s, t, tmp2));
			expectationvalues[i][i] = tmp.toFloat();
		}

		// clean basepairing matrices
		thissec = master.bottom;
		thissec.clearAllBp();
		while (thissec.next != null) {
			thissec = thissec.next;
			thissec.clearAllBp();
			thissec.clearAllExpectation();
		}

		ResultBundle result = new ResultBundle();
		result.structure = structure;
		result.reliability = reliability;
		result.basepairprob = basepairprob;
		result.singlebaseprob = singlebaseprob;
		result.expectationvalues = expectationvalues;
		result.entropyVal = entropyVal;
		result.entropyPercOfMax = entropyPercOfMax;
		result.entropyMax = entropyMax;
		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

//		System.out.print(nrdivisions + " - ");
		System.out.println("TOTAL TIME ELAPSED IN SCFG PART: "
				+ (int)((System.nanoTime() - starttime) * 1e-9) + " seconds ");
	//	System.out.println("            ...of which distributed: "
	//			+ (gridtime1 + gridtime2 + gridtime3) + " seconds ");
		//System.out.println();

		act.setProgress(1.0);
		return result;
	}

	private static double[][] createPhyloProb(Progress act, int userjobsnr,
			Tree tree, List<char[]> columns_char, List<String> names,
			final int length, Parameters param,
			AsynchronousJobExecutor executor, final boolean verbose, int execnr)
			throws InterruptedException {
		final long starttime = System.nanoTime();
		if(verbose){
		System.out.println("Timer (phylogeny) started. (time: "
				+ (System.nanoTime() - starttime) * 1e-9 + " s)");
		}
		if (verbose) {
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if (verbose) {
			System.out.println("Processing input...");
		}
		if (verbose) {
			act.setCurrentActivity("Evolutionary model: processing input");
		}
		List<int[]> columns = new ArrayList<int[]>();
		for (int i = 0; i < columns_char.size(); i++) {
			columns.add(MatrixTools.convertColumn(columns_char.get(i)));
		}
		if (verbose) {
			System.out.println("User wish for number of divisions: "
					+ userjobsnr);
		}

		// correct user input
		if (userjobsnr > length) {
			userjobsnr = length;
		} else if (userjobsnr < 1) {
			userjobsnr = 1;
		}
		int nrjobs = userjobsnr;

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if (verbose) {
			System.out
					.println("Actual nr. of divisions in phylogenetic calculations: "
							+ nrjobs);
		}
		if (verbose) {
			System.out.println("Generating matrices for each node...");
		}

		tree.generateLeafList(names); // to speed up calculations later

		final double[][] probmatrix = new double[length][length];
		List<PhyloJob> jobs = new ArrayList<PhyloJob>();
		// First generate "exp(Rt)" matrix for each node.
		// This is the same for all columns of that node.
		tree.getRoot().calculateChildrenMatrix(param.getSD(), param.getSV(),
				param.getSV1());

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out.println("Generating jobs...");
		}
		if (verbose) {
			act.setCurrentActivity("Evolutionary model: dividing tasks");
		}

		// Create jobs for SINGLE columns
		// now all single columns are in one job
		// Add the only single-column job
		int colcnt = 0;
		PhyloJob lastjob = new PhyloJob();
		lastjob.tree = Tree.copyTree(tree); // must copy the entire tree into
		// each phylojob
		lastjob.names = names; // must copy the names into each phylojob (for
		// debugging only)
		lastjob.startcol = 0;
		lastjob.jobid = 0;
		lastjob.endcol = length - 1;
		lastjob.type = false;
		lastjob.param = param;
		for (int col = 0; col < length; col++) {
			// send columns
			lastjob.columns.add(columns.get(col));
			colcnt++;
		}
		jobs.add(lastjob);

		// Create jobs for column PAIRS
		// count total number of column pairs
		int paircnt = 0;
		for (int i = 0; i < length; i++) {
			for (int j = i + 1; j < length; j++) {
				paircnt++;
			}
		}
		if(verbose){
			System.out.println("Total number of pairs: " + paircnt);
			System.out.println("Pairs in a job: " + (paircnt) / nrjobs);
		}
		// have to calculate matrices, but only once
		tree.getRoot().calculateChildrenMatrix(param.getDD(), param.getDV(),
				param.getDV1());

		int lastendcol = -1;
		int paircolcnt = 0;
		// add all jobs except last
		for (int jobnr = 0; jobnr < nrjobs - 1; jobnr++) {
			PhyloJob job = new PhyloJob();
			job.tree = Tree.copyTree(tree); // must copy the entire tree into
			// each phylojob
			job.names = names; // must copy the names into each phylojob
			job.type = true;
			job.param = param;
			job.jobid = jobnr + 1;

			int col1 = lastendcol + 1;
			job.startcol = col1;
			paircolcnt += length - col1;
			col1--;
			while (paircolcnt < ((long)(jobnr + 1) * (paircnt)) / nrjobs) {
				col1++;
				if (col1 == length) { // no more columns left, finish
					col1--;
					break;
				}
				job.columns.add(columns.get(col1));
				paircolcnt += length - col1;
			}

			job.endcol = col1;
			lastendcol = col1;
			for (int col2 = job.startcol + 1; col2 < length; col2++) {
				// send pairing columns
				job.columns2.add(columns.get(col2));
			}
			jobs.add(job);
			paircolcnt -= length - col1;
		}
		// Add the last job
		lastjob = new PhyloJob();
		lastjob.tree = Tree.copyTree(tree); // must copy the entire tree into
		// each phylojob
		lastjob.names = names; // must copy the names into each phylojob
		lastjob.startcol = lastendcol + 1;
		lastjob.endcol = length;
		lastjob.type = true;
		lastjob.param = param;
		lastjob.jobid = nrjobs;
		for (int col1 = lastjob.startcol; col1 < length; col1++) {
			// send columns
			lastjob.columns.add(columns.get(col1));
		}
		for (int col2 = lastjob.startcol + 1; col2 < length; col2++) {
			// send pairing columns
			lastjob.columns2.add(columns.get(col2));
		}
		jobs.add(lastjob);

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out
					.println("Total number of jobs in phylogenetic calculations: "
							+ jobs.size());
		}
		if (verbose) {
			System.out.println("Executing jobs...");
		}

		// start executing the jobs
		// counts how many phylojobs are done		
		final AtomicInteger finishedphylojobscount = new AtomicInteger(0); 

		act.setCurrentActivity("Evolutionary model: " +
					"calculating column probabilities");
		
		long gridstarttime = System.nanoTime();
		// execute single column jobs
		Progress singleColAct = act.getChildProgress(0.1);
		for (int jobnr = 0; jobnr < 1; jobnr++) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			PhyloJob job = jobs.get(jobnr);
			//final int jn = jobnr;
			final int startcol = job.startcol;

			final Progress jobAct = singleColAct
					.getChildProgress((float) job.columns.size()
							/ (float) colcnt);
			executor.startExecution(job, new JobListener() {
				public void jobFinished(JobResults result) {
				} // doesn't happen here

				public void jobFinished(List<ResultBundle> result) {
				} // doesn't happen here

				public void jobFinished(double[][] result) {
					for (int col = 0; col < result.length; col++) {
						probmatrix[col + startcol][col + startcol] = result[col][0];
					}
					finishedphylojobscount.incrementAndGet();
					jobAct.setProgress(1.0);
				}
			});
		}

		Progress doubleColAct = act.getChildProgress(0.9);

		// execute double column jobs
		for (int jobnr = 1; jobnr < jobs.size(); jobnr++) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			PhyloJob job = jobs.get(jobnr);
			//final int jn = jobnr;
			final int startcol = job.startcol;
			final int endcol = job.endcol;
			final Progress jobAct = doubleColAct
					.getChildProgress((float) job.columns.size()
							/ (float) length);
			if (job.columns2.size() != 0) {
				executor.startExecution(job, new JobListener() {
					public void jobFinished(JobResults result) {
					} // doesn't happen here

					public void jobFinished(List<ResultBundle> result) {
					} // doesn't happen here

					public void jobFinished(double[][] result) {
						for (int col1 = startcol; col1 < endcol + 1; col1++) {
							for (int col2 = col1 + 1; col2 < length; col2++) {
								probmatrix[col1][col2] = result[col1 - startcol][col2
										- startcol - 1];
								probmatrix[col2][col1] = probmatrix[col1][col2]; // make
								// it
								// symmetric
							}
						}
						finishedphylojobscount.incrementAndGet();
						jobAct.setProgress(1.0);
					}
				});
			} else {
				finishedphylojobscount.incrementAndGet();
				jobAct.setProgress(1.0);
			}
		}

		// wait for last thread to finish
		// wait for jobs to finish
		while (finishedphylojobscount.get() < jobs.size()) {
			Thread.sleep(100);
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
		}
		singleColAct.setProgress(1.0);
		doubleColAct.setProgress(1.0);
		act.setProgress(1.0);
		//act.setCurrentActivity("Evolutionary model: done");

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		System.out.println("TOTAL TIME ELAPSED IN PHYLOGENETIC PART: "
				+ (int)((System.nanoTime() - starttime) * 1e-9) + " seconds ");
		
		if(verbose){
			System.out.println("                    ...of which distributed: "
				+ (System.nanoTime() - gridstarttime) * 1e-9 + " seconds");
		}

		//System.out.println();

		// Result contains the a priori probability distribution matrix.
		return probmatrix;
	}

	static Sector findSector(int i, int j, Sector bottom) {
		// finds the sector in which the point i,j is located
		Sector current = bottom;
		// first find vertical position
		while (current.next != null) {
			if (i >= current.pos[0]
					&& (i < current.next.pos[0] || current.next.pos[0] == 0)) {
				break;
			}
			current = current.next;
		}
		// current now has the right column
		while (current.above != null) {
			if (j >= current.pos[1] - (i - current.pos[0])
					&& j < current.above.pos[1] - (i - current.pos[0])) {
				break;
			}
			current = current.above;
		}
		return current;

	}

	static int[] findPointST(int i, int j, Sector sector, int distance) {
		int result[] = new int[2];
		// lower half of the triangle, just return whatever sector is the one
		result[0] = sector.pos[0] + distance - i - 1;
		result[1] = j + i - (sector.pos[1] + sector.pos[0]);
		return result;
	}

	static int pairs(int a, int b) {
		int answer = 0;
		for (int c1 = 0; c1 < a; c1++) {
			for (int c2 = c1 + 1; c2 < b; c2++) {
				answer++;
			}
		}
		return answer;
	}
	
	static double log2(double val){
		return Math.log(val)/LOG_TWO;
	}

}