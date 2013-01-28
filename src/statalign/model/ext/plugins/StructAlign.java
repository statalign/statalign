package statalign.model.ext.plugins;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JToggleButton;

import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.MathArrays;

import statalign.base.InputData;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.hmm.Hmm;
import statalign.io.DataType;
import statalign.io.ProteinSkeletons;
import statalign.model.ext.ModelExtension;
import statalign.model.ext.McmcMove;
import statalign.model.ext.McmcCombinationMove;
import statalign.model.ext.GammaPrior;
import statalign.model.ext.plugins.structalign.*;
import statalign.model.ext.ParameterInterface;

import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.PluginParameters;

public class StructAlign extends ModelExtension implements ActionListener {
	
	/** The command line identifier of this plugin */
	//private static final String CMD_LINE_PLUGIN_ID = "structal";
	private final String pluginID = "structal";
	
	@Override
	public String getPluginID() {
		return pluginID;
	}
	
	JToggleButton myButton;
	
	public boolean globalSigma = true;
	double structTemp = 1;

	
	/** Alpha-C atomic coordinate for each sequence and each residue */
	public double[][][] coords;
	
	/** Alpha-C atomic coordinates under the current set of rotations/translations */
	public double[][][] rotCoords;
	
	/** Axis of rotation for each sequence */
	public double[][] axes;
	/** Rotation angle for each protein along the rotation axis */
	public double[] angles;
	/** Translation vector for each protein */
	public double[][] xlats;

	/** Parameters of structural drift */
	public double[] sigma2;
	public double sigma2Hier;
	public double nu;
	public double tau;
	public double epsilon;
	// TODO Allow starting values to be specified at command line/GUI
	
	/** Covariance matrix implied by current tree topology */
	public double[][] fullCovar;
	/** Current alignment between all leaf sequences */
	public String[] curAlign;
	
	/** Current log-likelihood contribution */
	public double curLogLike = 0;
	
	/** independence rotation proposal distribution */
	public RotationProposal rotProp;

	public double[][] oldCovar;
	public String[] oldAlign;
	public double oldLogLi;
	
	
	// TODO change the above public variables to package visible and put 
	// StructAlign.java in statalign.model.ext.plugins.structalign
	
	
	/** Priors */
	private double sigma2PriorShape = 0.001;
	private double sigma2PriorRate = 0.001;
	public GammaPrior sigma2Prior = new GammaPrior(sigma2PriorShape,sigma2PriorRate);
	
	private double tauPriorShape = 0.001;
	private double tauPriorRate = 0.001;
	public GammaPrior tauPrior = new GammaPrior(tauPriorShape,tauPriorRate);
	
	private double epsilonPriorShape = 0.001;
	private double epsilonPriorRate = 0.001;
	public GammaPrior epsilonPrior = new GammaPrior(epsilonPriorShape,epsilonPriorRate);
	
	private double sigma2HPriorShape = 0.001;
	private double sigma2HPriorRate = 0.001;
	public GammaPrior sigma2HPrior = new GammaPrior(sigma2HPriorShape,sigma2HPriorRate);
	
	private double nuPriorShape = 0.001;
	private double nuPriorRate = 0.001;
	public GammaPrior nuPrior = new GammaPrior(nuPriorShape,nuPriorRate);
	
	// priors for rotation and translation are uniform
	// so do not need to be included in M-H ratio
	
	/** Constant weights for 
	 * rotation/translation, sigma2, tau, sigma2Hier, nu, epsilon, subtree rotation, fixed-subtree alignment, and subtree rot+align combined */
	int[] paramPropWConst = { 0, 0, 3, 3, 3, 3, 2, 2, 2 };
	/** Weights per sequence for 
	 * rotation/translation, sigma2, tau, sigma2Hier, nu, epsilon, subtree rotation, fixed-subtree alignment, and subtree rot+align combined */
	int[] paramPropWPerSeq = { 5, 3, 0, 0, 0, 0, 0, 0, 0 };
	
	/** Total weights calculated as const+perseq*nseq */
	int[] paramPropWeights;
	/** Weights for proposing rotation vs translation vs library */
	int[] rotXlatWeights= { 25, 25, 10 };
	int[] subtreeRotXlatWeights = { 25, 25, 0};
//	int[] rotXlatWeights= { 25, 25, 0 };	// library off
	
	RotationMove rotationMove;
	
	
	/** Default proposal weights in this order: 
	 *  align, topology, edge, indel param, subst param, modelext param 
	*   { 35, 20, 15, 15, 10, 0 };
	*/
	private final int pluginProposalWeight = 50; 
	
	/** Proposal tuning parameters */
	public static final double angleP = 1000;
	public static final double xlatP = .1;
	
	public final double MIN_EPSILON = 0.1;
	
	@Override
	public List<JComponent> getToolBarItems() {
		myButton = new JToggleButton(new ImageIcon(ClassLoader.getSystemResource("icons/protein.png")));
    	myButton.setToolTipText("Structural alignment mode (for proteins only)");
    	myButton.addActionListener(this);
    	myButton.setEnabled(true);
    	myButton.setSelected(false);
    	return Arrays.asList((JComponent)myButton);
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		setActive(myButton.isSelected());
	}
	

	@Override
	public String getUsageInfo() {
		StringBuilder usage = new StringBuilder();
		usage.append("StructAlign version 1.0\n\n");
		usage.append("java -jar statalign.jar -plugin:structal=[OPTIONS]\n");
		
		return usage.toString();
	}

	@Override
	public void setActive(boolean active) {
		super.setActive(active);
		System.out.println("StructAlign plugin is now "+(active?"enabled":"disabled"));
	}
	
	@Override
	public void init(PluginParameters params) {
		if(params != null && params.getParameter(pluginID) != null) {
			// TODO parse plugin parameters
			setActive(true);
		}
	}
	
	@Override
	public void initRun(InputData inputData) throws IllegalArgumentException {
		HashMap<String, Integer> seqMap = new HashMap<String, Integer>();
		int i = 0;
		for(String name : inputData.seqs.seqNames)
			seqMap.put(name.toUpperCase(), i++);
		coords = new double[inputData.seqs.seqNames.size()][][];
		for(DataType data : inputData.auxData) {
			if(!(data instanceof ProteinSkeletons))
				continue;
			ProteinSkeletons ps = (ProteinSkeletons) data;
			for(i = 0; i < ps.names.size(); i++) {
				String name = ps.names.get(i).toUpperCase();
				if(!seqMap.containsKey(name))
					throw new IllegalArgumentException("structalign: missing sequence or duplicate structure for "+name);
				int ind = seqMap.get(name);
				int len = inputData.seqs.sequences.get(ind).replaceAll("-", "").length();
				List<double[]> cl = ps.coords.get(i);
				if(len != cl.size())
					throw new IllegalArgumentException("structalign: sequence length mismatch with structure file for seq "+name);
				coords[ind] = new double[len][];
				for(int j = 0; j < len; j++)
					 coords[ind][j] = Utils.copyOf(cl.get(j));
				RealMatrix temp = new Array2DRowRealMatrix(coords[ind]);
				RealVector mean = Funcs.meanVector(temp);
				for(int j = 0; j < len; j++)
					 coords[ind][j]= temp.getRowVector(j).subtract(mean).toArray();
				seqMap.remove(name);
			}
		}
		if(seqMap.size() > 0)
			throw new IllegalArgumentException("structalign: missing structure for sequence "+seqMap.keySet().iterator().next());
		
		rotProp = new RotationProposal(this);
		rotCoords = new double[coords.length][][];
		axes = new double[coords.length][];
		angles = new double[coords.length];
		xlats = new double[coords.length][];

		axes[0] = new double[] { 1, 0, 0 };
		angles[0] = 0;
		xlats[0] = new double[] { 0, 0, 0 };
						
		sigma2Hier = 1;
		nu = 1;
		tau = 5;
		epsilon = 5;
		
		// alternative initializations
		// actual initialization now occurs in beforeSampling()
		/*
		for(i = 1; i < axes.length; i++) {
			Transformation initial = rotProp.propose(i);
			axes[i] = initial.axis.toArray();
			angles[i] = initial.rot;
			xlats[i] = initial.xlat.toArray();
		}
		for(i = 0; i < axes.length; i++) {
			axes[i] = new double[] { 1, 0, 0 };
			angles[i] = 0;
			xlats[i] = new double[] { 0, 0, 0 };
		} */
		
		// number of branches in the tree is 2*leaves - 1
		if (globalSigma) {
			sigma2 = new double[1];
			paramPropWConst[3] = 0;
			paramPropWConst[4] = 0;
		}
		else {
			sigma2 = new double[2*coords.length - 1];
		}
		
		for(i = 0; i < sigma2.length; i++)
			sigma2[i] = 1;
		
		paramPropWeights = Utils.copyOf(paramPropWConst);
		for(i = 0; i < paramPropWeights.length; i++)
			paramPropWeights[i] += paramPropWPerSeq[i]*coords.length;
		
		/** Add alignment and rotation/translation moves */
		RotationMove rotationMove = new RotationMove(this,"rotation"); 
		addMcmcMove(rotationMove,10); // change to correct weight
		TranslationMove translationMove = new TranslationMove(this,"translation");
		addMcmcMove(translationMove,10); // change to correct weight
		LibraryMove libraryMove = new LibraryMove(this,"library");
		addMcmcMove(libraryMove,10); // change to correct weight
		AlignmentMove alignmentMove = new AlignmentMove(this,"alignment");
		addMcmcMove(alignmentMove,10); // change to correct weight
		
		ArrayList<McmcMove> alignmentRotation = new ArrayList<McmcMove>();
		alignmentRotation.add(alignmentMove);
		alignmentRotation.add(rotationMove);
		McmcCombinationMove alignmentRotationMove = 
			new McmcCombinationMove(alignmentRotation);
		addMcmcMove(alignmentRotationMove,10); // change to correct weight
		
		ArrayList<McmcMove> alignmentTranslation = new ArrayList<McmcMove>(); 
		alignmentTranslation.add(alignmentMove);
		alignmentTranslation.add(translationMove);
		McmcCombinationMove alignmentTranslationMove = 
			new McmcCombinationMove(alignmentTranslation);
		addMcmcMove(alignmentTranslationMove,10); // change to correct weight
		
		ArrayList<McmcMove> alignmentLibrary = new ArrayList<McmcMove>();
		alignmentLibrary.add(alignmentMove);
		alignmentLibrary.add(libraryMove);
		McmcCombinationMove alignmentLibraryMove = 
			new McmcCombinationMove(alignmentLibrary);
		addMcmcMove(alignmentLibraryMove,10); // change to correct weight
		
		/** Add moves for scalar parameters */
		StructAlignParameterInterface paramInterfaceGenerator = new StructAlignParameterInterface(this); 
		
		ParameterInterface tauInterface = paramInterfaceGenerator.new TauInterface();
		ContinuousPositiveParameterMove tauMove = 
			new ContinuousPositiveParameterMove(this,tauInterface,tauPrior,"tau");
		addMcmcMove(tauMove,10); // change to correct weight
		
		ParameterInterface epsilonInterface = paramInterfaceGenerator.new EpsilonInterface();
		ContinuousPositiveParameterMove epsilonMove = 
			new ContinuousPositiveParameterMove(this,epsilonInterface,epsilonPrior,"epsilon");
		epsilonMove.setMinValue(MIN_EPSILON);
		addMcmcMove(epsilonMove,10); // change to correct weight
				
		HierarchicalContinuousPositiveParameterMove sigma2HMove = null;
		HierarchicalContinuousPositiveParameterMove nuMove = null;
		if (!globalSigma) {
			ParameterInterface sigma2HInterface = paramInterfaceGenerator.new Sigma2HInterface();
			sigma2HMove = new HierarchicalContinuousPositiveParameterMove(this,sigma2HInterface,sigma2HPrior,"sigma2Hier");
			addMcmcMove(sigma2HMove,10); // change to correct weight
			
			ParameterInterface nuInterface = paramInterfaceGenerator.new NuInterface();
			nuMove = new HierarchicalContinuousPositiveParameterMove(this,nuInterface,nuPrior,"nu");
			addMcmcMove(nuMove,10); // change to correct weight
		}
		
		for (int j=0; j<sigma2.length; j++) {
			ParameterInterface sigma2Interface = paramInterfaceGenerator.new Sigma2Interface(j);
			ContinuousPositiveParameterMove m = new ContinuousPositiveParameterMove(
														this,sigma2Interface,
														sigma2Prior,"sigma2_"+i);
			addMcmcMove(m,10); // change to correct weight
			if (!globalSigma) {
				sigma2HMove.addChildMove(m);
				nuMove.addChildMove(m);
			}
		}
	}
	
	@Override
	public void beforeSampling(Tree tree) {
		Funcs.initLSRotations(tree,coords,xlats,axes,angles);
	}
	
	
	@Override
	public double logLikeFactor(Tree tree) {
		String[] align = tree.getState().getLeafAlign();
		checkConsAlign(align); 		
		curAlign = align;
		
		double[][] covar = calcFullCovar(tree);
		checkConsCovar(covar); 
		fullCovar = covar;
		
		if(!checkConsRots() && rotCoords[0] == null)
			calcAllRotations();
		
		/** TESTING
		
		System.out.println();
		System.out.println("Parameters for structural log likelihood:");		
		System.out.println("Sigma2: " + sigma2);
		System.out.println("Theta: " + theta);
		System.out.println("Branch length: " + (tree.root.left.edgeLength + tree.root.right.edgeLength));
		
		System.out.println("Rotation matrices:");
		for(int i = 1; i < xlats.length; i++) {
			Rotation rot = new Rotation(new Vector3D(axes[i]), angles[i]);
			double[][] m = rot.getMatrix();
			for(int j = 0; j < m.length; j++)
				System.out.println(Arrays.toString(m[j]));
		}
		
		System.out.println("Translations:");
		for(int i = 0; i < xlats.length; i++)
			System.out.println(Arrays.toString(xlats[i]));
		
		/** END TESTING */
		
		
		double logli = calcAllColumnContrib();
		checkConsLogLike(logli); 
		curLogLike = logli;
		
		// testing
		//System.out.println("Total log likelihood " + curLogLike);
		
		return curLogLike;
	}
	
	public double calcAllColumnContrib() {
		String[] align = curAlign;
		double logli = 0;
		int[] inds = new int[align.length];		// current char indices
		int[] col = new int[align.length];  
		for(int i = 0; i < align[0].length(); i++) {
			for(int j = 0; j < align.length; j++)
				col[j] = align[j].charAt(i) == '-' ? -1 : inds[j]++;
			double ll = columnContrib(col); 
			logli += ll;
			//System.out.println("Column: " + Arrays.toString(col) + "  ll: " + ll);
		}
		return structTemp * logli;
	}
	// TODO Change visibility of this to package, after moving
	// StructAlign.java to statalign.model.ext.plugins.structalign

	private boolean checkConsAlign(String[] align) {
		if(!Utils.DEBUG || curAlign == null)
			return false;
		if(align.length != curAlign.length)
			throw new Error("Inconsistency in StructAlign, alignment length: "+align.length+", "+curAlign.length);
		for(int i = 0; i < align.length; i++)
			if(!align[i].equals(curAlign[i]))
				throw new Error("Inconsistency in StructAlign, alignment: "+align[i]+", "+curAlign[i]);
		return true;
	}

	private boolean checkConsCovar(double[][] covar) {
		if(!Utils.DEBUG || fullCovar == null)
			return false;
		if(covar.length != fullCovar.length)
			throw new Error("Inconsistency in StructAlign, covar matrix length: "+covar.length+", "+fullCovar.length);
		for(int i = 0; i < covar.length; i++) {
			if(covar[i].length != fullCovar[i].length)
				throw new Error("Inconsistency in StructAlign, covar matrix "+i+" length: "+covar[i].length+", "+fullCovar[i].length);
			for(int j = 0; j < covar[i].length; j++)
				if(Math.abs(covar[i][j]-fullCovar[i][j]) > 1e-5)
					throw new Error("Inconsistency in StructAlign, covar matrix "+i+","+j+" value: "+covar[i][j]+", "+fullCovar[i][j]+", "+tau+", "+epsilon);
		}
		return true;
	}
	
	private boolean checkConsRots() {
		if(!Utils.DEBUG || rotCoords[0] == null)
			return false;
		double[][][] rots = new double[rotCoords.length][][];
		for(int i = 0; i < rots.length; i++) {
			rots[i] = new double[rotCoords[i].length][];
			for(int j = 0; j < rots[i].length; j++)
				rots[i][j] = MathArrays.copyOf(rotCoords[i][j]);
		}
		calcAllRotations();
		for(int i = 0; i < rots.length; i++)
			for(int j = 0; j < rots[i].length; j++)
				for(int k = 0; k < rots[i][j].length; k++)
					if(Math.abs(rots[i][j][k]-rotCoords[i][j][k]) > 1e-5)
						throw new Error("Inconsistency in StructAlign, rotation "+i+","+j+","+k+": "+rots[i][j][k]+" vs "+rotCoords[i][j][k]);
		return true;
	}

	private boolean checkConsLogLike(double logli) {
		if(!Utils.DEBUG || curLogLike == 0)
			return false;
		if(Math.abs(logli-curLogLike) > 1e-5)
			throw new Error("Inconsistency in StructAlign, log-likelihood "+logli+" vs "+curLogLike);
		return true;
	}

	/**
	 * Calculates the structural likelihood contribution of a single alignment column
	 * @param col the column, id of the residue for each sequence (or -1 if gapped in column)
	 * @return the likelihood contribution
	 */
	public double columnContrib(int[] col) {
		// count the number of ungapped positions in the column
		int numMatch = 0;
		for(int i = 0; i < col.length; i++){
			if(col[i] != -1)
				numMatch++;
		}
		if(numMatch == 0) 
			return 1;
		// collect indices of ungapped positions
		int[] notgap = new int[numMatch];
		int j = 0;
		for(int i = 0; i < col.length; i++)
			if(col[i] != -1)
				notgap[j++] = i;
		
		// extract covariance corresponding to ungapped positions
		double[][] subCovar = Funcs.getSubMatrix(fullCovar, notgap, notgap);
		// create normal distribution with mean 0 and covariance subCovar
		MultiNormCholesky multiNorm = new MultiNormCholesky(new double[numMatch], subCovar);
		
		double logli = 0;
		double[] vals = new double[numMatch];
		// loop over all 3 coordinates
		
		/*System.out.println("Calculating log likelihood: ");
		System.out.println("Mean: " + Arrays.toString(multiNorm.getMeans()));
		System.out.println("Variance: " + Arrays.toString(subCovar[0]));*/
		for(j = 0; j < 3; j++){
			for(int i = 0; i < numMatch; i++)
				vals[i] = rotCoords[notgap[i]][col[notgap[i]]][j];
			//System.out.println("Values: " + Arrays.toString(vals));
			logli += multiNorm.logDensity(vals);
			//System.out.println("LL: " + multiNorm.logDensity(vals));
		}
		return logli;
	}

	/**
	 * extracts the specified rows and columns of a 2d array
	 * @param matrix, 2d array from which to extract; rows, rows to extract; cols, columns to extract
	 * @return submatrix
	 */
		
	private void calcAllRotations() {
		for(int i = 0; i < coords.length; i++)
			calcRotation(i);
	}
	
	public void calcRotation(int ind) {
		double[][] ci = coords[ind], rci = rotCoords[ind];
		if(rci == null)
			rci = rotCoords[ind] = new double[ci.length][];
		Rotation rot = new Rotation(new Vector3D(axes[ind]), angles[ind]);
		for(int i = 0; i < ci.length; i++) {
			rci[i] = rot.applyTo(new Vector3D(ci[i])).add(new Vector3D(xlats[ind])).toArray();
		}
	}
	// TODO Change visibility of this to package, after moving
	// StructAlign.java into statalign.model.ext.plugins.structalign.

	/**
	 * return the full covariance matrix for the tree topology and branch lengths
	 */	
	public double[][] calcFullCovar(Tree tree) {
		// I'm assuming that tree.names.length is equal to the number of vertices here
		double[][] distMat = new double[tree.names.length][tree.names.length];
		calcDistanceMatrix(tree.root, distMat);
		//System.out.print("Distance: " + distMat[0][1]);
		
		//System.out.println("Current tree:");
		//printTree(tree.root, "o");
		
		for(int i = 0; i < tree.names.length; i++)
			for(int j = i; j < tree.names.length; j++)
				distMat[j][i] = distMat[i][j] = tau * Math.exp(-distMat[i][j]);
		for(int i = 0; i < tree.names.length; i++)
			distMat[i][i] += epsilon;
		return distMat;
	}
	

	public void printTree(Vertex v, String vname){
		System.out.println(vname +"-" + v.name + ": " + v.edgeLength);
		if(v.left!=null){
			printTree(v.left, vname + "l");
			printTree(v.right, vname + "r");
		}
	}
	
	
	/**
	 * recursive algorithm to traverse tree and calculate distance matrix between leaves 
	 */		
	public int[] calcDistanceMatrix(Vertex vertex, double[][] distMat){
		int[] subTree = new int[distMat.length + 1];
		
		// either both left and right are null or neither is
		if(vertex.left != null){
			int[] subLeft  = calcDistanceMatrix(vertex.left, distMat);
			int[] subRight = calcDistanceMatrix(vertex.right, distMat);
			int i = 0;
			while(subLeft[i] > -1){
				subTree[i] = subLeft[i];
				i++;
			}
			for(int j = 0; i+j < subTree.length; j++)
				subTree[i+j] = subRight[j];
		}
		else{
			subTree[0] = vertex.index;
			for(int j = 1; j < subTree.length; j++)
				subTree[j] = -1;
		}

		if (globalSigma) {
			addEdgeLength(distMat, subTree, vertex.edgeLength * sigma2[0] / (2*tau));	
		}
		else {
			addEdgeLength(distMat, subTree, vertex.edgeLength * sigma2[vertex.index] / (2*tau));
		}
		/*System.out.println();
		System.out.println("Distmat:");
		for(int i = 0; i < distMat.length; i++)
			for(int j = 0; j < distMat[0].length; j++)
				System.out.println(distMat[i][j]);*/
		return subTree;
	}
		
	// adds the length of the current edge to the distance between all leaves
	// of a subtree to all other leaves
	// 'rows' contains the indices of vertices in the subtree
	public void addEdgeLength(double[][] distMat, int[] subTree, double edgeLength){
		
		int i = 0;
		while(subTree[i] > -1){
			for(int j = 0; j < distMat.length; j++){  
				distMat[subTree[i]][j] += edgeLength;
				distMat[j][subTree[i]] += edgeLength;
			}
			i++;		
		}
			
		// edge length should not be added to distance between vertices in the subtree
		// subtract the value from these entries of the distance matrix
		i = 0;
		while(subTree[i] > -1){
			int j = 0;
			while(subTree[j] > -1){
				distMat[subTree[i]][subTree[j]] -= edgeLength;
				distMat[subTree[j]][subTree[i]] -= edgeLength;
				j++;
			}
			i++;
		}
	}


	@Override
	public int getParamChangeWeight() {
		// TODO test converge and tune value
		return pluginProposalWeight;
	}

	@Override
	public double logLikeModExtParamChange(Tree tree, ModelExtension ext) {
		// current log-likelihood always precomputed (regardless of whether ext == this)
		return curLogLike;
	}
	
	@Override
	public double logLikeAlignChange(Tree tree, Vertex selectRoot) {
		oldAlign = curAlign;
		oldLogLi = curLogLike;
		curAlign = tree.getState().getLeafAlign();
		curLogLike = calcAllColumnContrib();
		return curLogLike;
	}
	
	@Override
	public void afterAlignChange(Tree tree, Vertex selectRoot, boolean accepted) {
		if(accepted)	// accepted, do nothing
			return;
		// rejected, restore
		curAlign = oldAlign;
		curLogLike = oldLogLi;
	}
	
	@Override
	public double logLikeTreeChange(Tree tree, Vertex nephew) {
		oldCovar = fullCovar;
		oldAlign = curAlign;
		oldLogLi = curLogLike;
		fullCovar = calcFullCovar(tree);
		curAlign = tree.getState().getLeafAlign();
		curLogLike = calcAllColumnContrib();
		return curLogLike;
	}
	
	@Override
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {
		if(accepted)	// accepted, do nothing
			return;
		// rejected, restore
		fullCovar = oldCovar;
		curAlign = oldAlign;
		curLogLike = oldLogLi;
	}
	
	@Override
	public double logLikeEdgeLenChange(Tree tree, Vertex vertex) {
		// do exactly the same as for topology change
		return logLikeTreeChange(tree, vertex);
	}
	
	@Override
	public void afterEdgeLenChange(Tree tree, Vertex vertex, boolean accepted) {
		// do exactly the same as for topology change
		afterTreeChange(tree, vertex, accepted);
	}
	
	@Override
	public double logLikeIndelParamChange(Tree tree, Hmm hmm, int ind) {
		// does not affect log-likelihood
		return curLogLike;
	}
	
	@Override
	public double logLikeSubstParamChange(Tree tree, SubstitutionModel model,
			int ind) {
		// does not affect log-likelihood
		return curLogLike;
	}

	
	/** Rotation Library
	 * 
	 * @author Challis
	 *
	 */
	
	
	
	
	
	
	
	
	/**
	 * 
	 * @author Challis
	 *
	 */
	
	
	/**
	 * For storing helpful functions
	 * @author Challis
	 *
	 */
	
	public double[][] getRowSub(double[][] coord, ArrayList<Integer> rows){
		double[][] sub = new double[rows.size()][coord[0].length];
		
		for(int i = 0; i < sub.length; i++)
			for(int j = 0; j < sub[0].length; j++)
				sub[i][j] = coord[rows.get(i)][j];
		return sub;
	}

	// </StructAlign>
}


