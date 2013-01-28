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
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.MathArrays;

import statalign.base.InputData;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.hmm.Hmm;
import statalign.io.DataType;
import statalign.io.ProteinSkeletons;
import statalign.model.ext.ModelExtension;
import statalign.model.ext.plugins.structalign.*;

import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.PluginParameters;
import statalign.utils.BetaDistribution;
import statalign.utils.GammaDistribution;
import statalign.utils.NormalDistribution;
import cern.jet.math.Bessel;

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
	public double sigma2Hier = 1;
	public double nu = 1;
	public double tau = 5;
	public double epsilon = 5;
	// TODO Allow starting values to be specified at command line/GUI
	
	/** Covariance matrix implied by current tree topology */
	double[][] fullCovar;
	/** Current alignment between all leaf sequences */
	private String[] curAlign;
	
	/** Current log-likelihood contribution */
	public double curLogLike = 0;
	
	/** independence rotation proposal distribution */
	public RotationProposal rotProp;
	
	/** MCMC moves */
	addMcmcMove(new RotationMove(this));
	TranslationMove translationMove = new TranslationMove(this);
	LibraryMove libraryMove = new LibraryMove(this);
	
	
	// TODO change the above public variables to package visible and put 
	// StructAlign.java in statalign.model.ext.plugins.structalign
	
	private double[][] oldCovar;
	private String[] oldAlign;
	private double oldLogLi;
	
	
	
	
	/** Priors */
	// tau - gamma prior, uses shape/scale parameterization
	GammaDistribution tauPrior = new GammaDistribution(1, 100);
	
	// epsilon - gamma prior
	GammaDistribution epsilonPrior = new GammaDistribution(1, 100);
	
	// sigma2Hier - gamma prior
	GammaDistribution sigma2HPrior = new GammaDistribution(1, 100);
	
	// nu - gamma prior
	GammaDistribution nuPrior = new GammaDistribution(1, 100);
	
	// initialize hierarchical distribution for sigma2
	GammaDistribution sigma2HierDist = new GammaDistribution(nu, sigma2Hier / nu);

	// priors for rotation and translation are uniform
	// so do not need to be included in M-H ratio
	
	/** Constant weights for 
	 * rotation/translation, sigma2, tau, sigma2Hier, nu, epsilon, subtree rotation, fixed-subtree alignment, and subtree rot+align combined */
	int[] paramPropWConst = { 0, 0, 3, 3, 3, 3, 2, 2, 2 };
	/** Weights per sequence for 
	 * rotation/translation, sigma2, tau, sigma2Hier, nu, epsilon, subtree rotation, fixed-subtree alignment, and subtree rot+align combined */
	int[] paramPropWPerSeq = { 5, 3, 0, 0, 0, 0, 0, 0, 0 };
	
	private int[] mcmcMoveWeights;


	/** Total weights calculated as const+perseq*nseq */
	int[] paramPropWeights;
	/** Weights for proposing rotation vs translation vs library */
	int[] rotXlatWeights= { 25, 25, 10 };
	int[] subtreeRotXlatWeights = { 25, 25, 0};
//	int[] rotXlatWeights= { 25, 25, 0 };	// library off
	
	
	/** Default proposal weights in this order: 
	 *  align, topology, edge, indel param, subst param, modelext param 
	*   { 35, 20, 15, 15, 10, 0 };
	*/
	private final int pluginProposalWeight = 50; 
	
	/** Proposal tuning parameters */
	// higher values lead to smaller step sizes
	//proposalCounts = new int[6];
	//proposalCounts = new int[6];
	//proposalCounts = new int[6];
	
	//private static final double tauP = 10;
	//private static final double sigma2P = 10;
	//private static final double sigma2HP = 10;
	//private static final double nuP = 10;
	//private static final double axisP = 100;
	public static final double angleP = 1000;
	// higher values lead to bigger step sizes
	public static final double xlatP = .1;
	
	public double MIN_EPSILON = 2;
	
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
		
		int others = 4;
		// number of branches in the tree is 2*leaves - 1
		if (globalSigma) {
			sigma2 = new double[1];
			paramPropWConst[3] = 0;
			paramPropWConst[4] = 0;
		}
		else {
			sigma2 = new double[2*coords.length - 1];
		}
		
		if (globalSigma) {
			proposalCounts = new int[3]; // Includes tau & epsilon
			acceptanceCounts = new int[3];
			proposalWidthControlVariables = new double[3];
		}
		else {
			proposalCounts = new int[sigma2.length+others]; // Includes tau, sigma2Hier, nu, and epsilon
			acceptanceCounts = new int[sigma2.length+others];
			proposalWidthControlVariables = new double[sigma2.length+others];			
		}
		for (int ii=0; ii<proposalCounts.length; ii++) {
			proposalCounts[ii] = 0;
			acceptanceCounts[ii] = 0;
			proposalWidthControlVariables[ii] = 1; 
			// This needs to be a variable that, when bigger, increases the 
			// width of the proposal.
		}
		
		for(i = 0; i < sigma2.length; i++)
			sigma2[i] = 1;
		
		// Indices of the variables in the arrays of proposal counts and acceptance counts.
		tauInd = sigma2.length;
		epsilonInd = sigma2.length+1;
		sigma2HInd = sigma2.length+2;
		nuInd = sigma2.length+3;
		
		
		tau = 5;
		tauProposed = 0;
		tauAccept = 0;
		rotProposed = 0;
		rotAccept = 0;
		xlatProposed = 0;
		xlatAccept = 0;
		libProposed = 0;
		libAccept = 0;
		
		paramPropWeights = Utils.copyOf(paramPropWConst);
		for(i = 0; i < paramPropWeights.length; i++)
			paramPropWeights[i] += paramPropWPerSeq[i]*coords.length;

	}
	
	@Override
	public void beforeSampling(Tree tree) {
		initLSRotations(tree);
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
		double[][] subCovar = getSubMatrix(fullCovar, notgap, notgap);
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
	public double[][] getSubMatrix(double[][] matrix, int[] rows, int[] cols) {
		double[][] submat = new double[rows.length][cols.length];
		for(int i = 0; i < rows.length; i++)
			for(int j = 0; j < cols.length; j++)
				submat[i][j] = matrix[rows[i]][cols[j]];
		return submat;
	}
	
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
	public void proposeParamChange(Tree tree) {
		int param = Utils.weightedChoose(paramPropWeights);
		if(param == 0) {
			writeRotationFiles(tree);
			// proposing rotation/translation of a single sequence
			int rotxlat = Utils.weightedChoose(rotXlatWeights);
			if (rotxlat == 0) {
				rotationMove.move();
			}
			else if (rotxlat == 1) {
				translationMove.move();
			}
			else if (rotxlat == 2) {
				libraryMove.move();
			}
			
			// choose sequence to rotate/translate (never rotate 1st one)
			int omit = (rotxlat != 1 ? 1 : 0);
			int ind = Utils.generator.nextInt(coords.length - omit) + omit;
			
			double[] oldax = MathArrays.copyOf(axes[ind]);
			double oldang = angles[ind];
			double[] oldxlat = MathArrays.copyOf(xlats[ind]);
			double[][] oldrots = rotCoords[ind];
			rotCoords[ind] = null;	// so that calcRotation creates new array
			double oldll = curLogLike;

			double logProposalRatio = 0;

			switch(rotxlat) {
			case 0:
				// rotation of a single sequence
				rotProposed++;
				// axes[ind] = vonMisesFisher.simulate(axisP, new ArrayRealVector(axes[ind])).toArray();
				double smallAngle = vonMises.simulate(angleP, 0);
				
				RealVector randomAxis = new ArrayRealVector(3);
				for(int i = 0; i < 3; i++)
					randomAxis.setEntry(i, Utils.generator.nextGaussian());
				randomAxis.unitize();
				
				Rotation Q = new Rotation(new Vector3D(randomAxis.toArray()), smallAngle);
				Rotation R = new Rotation(new Vector3D(axes[ind]), angles[ind]);
				
				R = Q.applyTo(R);
			
				axes[ind] = R.getAxis().toArray();
				angles[ind] = R.getAngle();
				
				// logProposalRatio is 0 because prior is uniform and proposal is symmetric
				
				break;
			case 1:
				// translation of a single sequence
				xlatProposed++;
				for(int i = 0; i < 3; i++)
					xlats[ind][i] = Utils.generator.nextGaussian() * xlatP + xlats[ind][i];  
				
				// logProposalRatio is 0 because prior is uniform and proposal is symmetric
				
				break;
			case 2:
				// library proposal of a single sequence
				libProposed++;
				Transformation old = new Transformation(axes[ind], angles[ind], xlats[ind]);
				// transformation should be relative to reference protein
				old.xlat = old.xlat.subtract(new ArrayRealVector(xlats[0]));
				Transformation prop = rotProp.propose(ind);
				axes[ind] = prop.rotMatrix.getAxis().toArray();
				angles[ind] = prop.rotMatrix.getAngle();
				xlats[ind] = prop.xlat.toArray();

				// library density 
				logProposalRatio = rotProp.libraryLogDensity(ind, old) - 
						  rotProp.libraryLogDensity(ind, prop);
				
				// proposed translation is relative to reference protein
				for(int i = 0; i < 3; i++)
					xlats[ind][i] += xlats[0][i];
							
				break;
			}

			calcRotation(ind);
			curLogLike = calcAllColumnContrib();
			if(isParamChangeAccepted(logProposalRatio)) {
				// accepted, nothing to do
				switch(rotxlat) {
				case 0:
					rotAccept++;
					break;
				case 1:
					xlatAccept++;
					break;
				case 2:
					libAccept++;
					break;
				}
//				if(Utils.DEBUG)
//					System.out.println(new String[] { "rot", "xlat", "library" }[rotxlat]+" accepted");
			} else {
				// rejected, restore
//				if(Utils.DEBUG)
//					System.out.println(new String[] { "rot", "xlat", "library" }[rotxlat]+" rejected");
				axes[ind] = oldax;
				angles[ind] = oldang;
				xlats[ind] = oldxlat;
				rotCoords[ind] = oldrots;
				curLogLike = oldll;
			}
				
		} else if(param == 1 || param == 2){
			
			int sigmaInd = 0;
			if(param == 1 && !globalSigma) {	// select sigma to propose if not global
				sigmaInd = Utils.generator.nextInt(2*coords.length-2);
				if(sigmaInd >= tree.root.index)
					sigmaInd++;
			}
			
			// proposing new sigma
			double oldpar = param == 1 ? sigma2[sigmaInd] : tau;
			double[][] oldcovar = fullCovar;
			double oldll = curLogLike;
			
			double logProposalRatio = 0;
			
			GammaDistribution proposal;
			GammaDistribution reverse;
			
			if(param == 1){
				//System.out.println("sigma2["+sigmaInd+"] = "+sigma2[sigmaInd]);
				//sigProposed[sigmaInd]++;
				proposalCounts[sigmaInd]++;
				double sigma2P = 1.0d/proposalWidthControlVariables[sigmaInd];
				proposal = new GammaDistribution(sigma2P + 0.001, oldpar / sigma2P + 0.001);
				//proposal = new GammaDistribution((0.001+oldpar)/sigma2P, sigma2P);
				sigma2[sigmaInd] = proposal.sample();
				reverse = new GammaDistribution(sigma2P + 0.001, sigma2[sigmaInd] / sigma2P + 0.001);
				//reverse = new GammaDistribution((0.001+sigma2[sigmaInd])/sigma2P, sigma2P);				
			} else{
				//tauProposed++;
				proposalCounts[tauInd]++;
				// creates a gamma distribution with mean theta & variance controlled by thetaP
				double tauP = 1.0d/proposalWidthControlVariables[tauInd];
				proposal = new GammaDistribution(tauP + 0.001, oldpar / tauP + 0.001);
				//proposal = new GammaDistribution((0.001+oldpar)/tauP, tauP);
				tau = proposal.sample();
				reverse = new GammaDistribution(tauP + 0.001, tau / tauP + 0.001);
				//reverse = new GammaDistribution((0.001+tau)/tauP, tauP);
			}
			// TODO do not recalculate distances, only the covariance matrix
			fullCovar = calcFullCovar(tree);
			curLogLike = calcAllColumnContrib();
			
			if(param == 1){
				logProposalRatio = Math.log(reverse.density(oldpar)) 
				            - Math.log(proposal.density(sigma2[sigmaInd]));
				if (!globalSigma) {
					logProposalRatio += Math.log(sigma2HierDist.density(sigma2[sigmaInd]))
							    - Math.log(sigma2HierDist.density(oldpar));
				}
			}
			else {
				logProposalRatio = Math.log(tauPrior.density(tau)) + Math.log(reverse.density(oldpar)) 
				- Math.log(tauPrior.density(oldpar)) - Math.log(proposal.density(tau));
//				System.out.println("Tau proposed");
//				System.out.println("oldll: " + oldll);
//				System.out.println("curLogLike: " + curLogLike);
//				System.out.println("logProposalRatio: " + logProposalRatio);
				
			}
			if (isParamChangeAccepted(logProposalRatio)) {
				if (param == 1)
					//sigAccept[sigmaInd]++;
					acceptanceCounts[sigmaInd]++;
				else {
					//tauAccept++;
					acceptanceCounts[tauInd]++;
					// accepted, nothing to do
					//System.out.println("tau accepted\n");
				}
			} else {
				// rejected, restore
				if(param == 1)
					sigma2[sigmaInd] = oldpar;
				else{
					//System.out.println("tau rejected\n");
					tau = oldpar;
				}
					
				fullCovar = oldcovar;
				curLogLike = oldll;
			}
		} else if(param == 3 || param == 4){

			// proposing new sigma2Hier/nu
			double oldpar = param == 3 ? sigma2Hier : nu;
			double oldsigll = 0;
			for(int i = 0; i < tree.vertex.length-1; i++)
				if(i != tree.root.index)
					oldsigll += Math.log(sigma2HierDist.density(sigma2[i]));
			
			double logProposalRatio = 0;
			
			GammaDistribution proposal;
			GammaDistribution reverse;
			
			if(param == 3){
				//sigHProposed++;
				proposalCounts[sigma2HInd]++;
				double sigma2HP = 1.0d/proposalWidthControlVariables[sigma2HInd];
				proposal = new GammaDistribution(sigma2HP + 0.001, oldpar / sigma2HP + 0.001);
				//proposal = new GammaDistribution((0.001+oldpar)/sigma2HP, sigma2HP);
				sigma2Hier = proposal.sample();
				reverse = new GammaDistribution(sigma2HP + 0.001, sigma2Hier / sigma2HP + 0.001);
				//reverse = new GammaDistribution((0.001+sigma2Hier)/sigma2HP, sigma2HP);
			} else{
				//nuProposed++;
				proposalCounts[nuInd]++;
				double nuP = 1.0d/proposalWidthControlVariables[nuInd];
				// creates a gamma distribution with mean nu & variance controlled by nuP
				proposal = new GammaDistribution(nuP + 0.001, oldpar / nuP + 0.001);
				//proposal = new GammaDistribution((0.001+oldpar)/nuP, nuP);
				nu = proposal.sample();
				reverse = new GammaDistribution(nuP + 0.001, nu / nuP + 0.001);
				//reverse = new GammaDistribution((0.001+nu)/nuP, nuP);
			}
			
			sigma2HierDist = new GammaDistribution(nu, sigma2Hier / nu);
			
			double newsigll = 0;
			for(int i = 0; i < tree.vertex.length-1; i++)
				if(i != tree.root.index)
					newsigll += Math.log(sigma2HierDist.density(sigma2[i]));
					
			if(param == 3)
				logProposalRatio = newsigll + Math.log(sigma2HPrior.density(sigma2Hier)) + Math.log(reverse.density(oldpar)) 
				- oldsigll - Math.log(sigma2HPrior.density(oldpar)) - Math.log(proposal.density(sigma2Hier));
			else
				logProposalRatio = newsigll + Math.log(nuPrior.density(nu)) + Math.log(reverse.density(oldpar)) 
				- oldsigll - Math.log(nuPrior.density(oldpar)) - Math.log(proposal.density(nu));
			
			if(isParamChangeAccepted(logProposalRatio)) {
				if(param == 3) {
					sigHAccept++;
					acceptanceCounts[sigma2HInd]++;
				}
				else {
					nuAccept++;
					acceptanceCounts[nuInd]++;
				}
				// accepted, nothing to do
			} else {
				// rejected, restore
				if(param == 3)
					sigma2Hier = oldpar;
				else
					nu = oldpar;
			}
		} else if(param == 5){
			
			// proposing new epsilon
			double oldpar = epsilon;
			double[][] oldcovar = fullCovar;
			double oldll = curLogLike;
			
			double logProposalRatio = 0;
			
			GammaDistribution proposal;
			GammaDistribution reverse;
			
			proposalCounts[epsilonInd]++;
			double epsilonP = 1.0d/proposalWidthControlVariables[epsilonInd];
			proposal = new GammaDistribution(epsilonP + .001, oldpar / epsilonP+.001);
			epsilon = proposal.sample();
			reverse = new GammaDistribution(epsilonP +.001, epsilon / epsilonP+.001);
			
			// TODO do not recalculate distances, only the covariance matrix
			fullCovar = calcFullCovar(tree);
			// for(int i = 0; i < coords.length; i++)
			//	fullCovar[i][i] = fullCovar[i][i] - oldpar + epsilon;
			curLogLike = calcAllColumnContrib();
			
			if(param == 1)
				logProposalRatio =  + Math.log(reverse.density(oldpar)) 
				            - Math.log(epsilonPrior.density(oldpar)) - Math.log(proposal.density(epsilon));

			if(epsilon > MIN_EPSILON && isParamChangeAccepted(logProposalRatio)) {
					acceptanceCounts[epsilonInd]++;
				// accepted, nothing to do
			} else {
				// rejected, restore
				epsilon = oldpar;
				fullCovar = oldcovar;
				curLogLike = oldll;
			}
			
		} else if(param==6 || param==7 || param==8) { // propose subtree move
			
			boolean proposingFixedSubtreeRotation = false;
			boolean proposingFixedSubtreeAlignment = false;
			System.out.print("Subtree ");
			if (param==6) {
				proposingFixedSubtreeRotation = true;
				System.out.print("rotation ");
				subtreeRotProposed++;
			}
			else if (param==7){
				proposingFixedSubtreeAlignment = true;
				System.out.print("alignment ");
				subtreeAlignProposed++;
			}
			else {
				proposingFixedSubtreeRotation = true;
				proposingFixedSubtreeAlignment = true;
				System.out.print("rotation + alignment");
				subtreeRotAlignProposed++;
			}
			
			double[][] oldaxes = null;
			double[] oldangles = null;
			double[][] oldxlats = null;
			double[][][] oldrots = null;
			
			double oldll = curLogLike;
			double logProposalRatio = 0;
			
			Vertex subtreeRoot = sampleVertex(tree);
			ArrayList<Integer> subtreeLeaves = collectLeaves(subtreeRoot);
			if(subtreeLeaves.contains(0)){	// check if subtree contains reference protein
				ArrayList<Integer> complement = new ArrayList<Integer>(0);
				for(int i = 0; i < coords.length; i++)
					complement.add(i);
				for(int i = 0; i < subtreeLeaves.size(); i++) {
					complement.remove(subtreeLeaves.get(i));
				}
				subtreeLeaves = complement;
			}
			System.out.print(" {");
			for(int i = 0; i < subtreeLeaves.size(); i++){
				int j = subtreeLeaves.get(i);
				System.out.print(j);
				if (i < (subtreeLeaves.size() - 1)) {
					System.out.print(",");
				}
			}
			System.out.print("} ");

			int index = subtreeLeaves.get(Utils.generator.nextInt(subtreeLeaves.size()));


			if (proposingFixedSubtreeRotation) {
				oldaxes = new double[axes.length][axes[0].length];
				oldangles = new double[angles.length];
				oldxlats = new double[xlats.length][xlats.length];
				oldrots = new double[rotCoords.length][rotCoords[0].length][rotCoords[0][0].length];

				for(int i = 0; i < subtreeLeaves.size(); i++){
					int j = subtreeLeaves.get(i);
					oldaxes[j] = MathArrays.copyOf(axes[j]);
					oldangles[j] = angles[j];
					oldxlats[j] = MathArrays.copyOf(xlats[j]);
					oldrots[j] = rotCoords[j];
					rotCoords[j] = null;	// so that calcRotation creates new array
				}
				
				int rotxlat = Utils.weightedChoose(subtreeRotXlatWeights);
				
				switch(rotxlat) {
				case 0:
					System.out.print("Rotation: ");
					// rotation of group of proteins
					rotProposed++;
					// axes[ind] = vonMisesFisher.simulate(axisP, new ArrayRealVector(axes[ind])).toArray();
					double smallAngle = vonMises.simulate(angleP, 0);
					
					RealVector randomAxis = new ArrayRealVector(3);
					for(int i = 0; i < 3; i++)
						randomAxis.setEntry(i, Utils.generator.nextGaussian());
					randomAxis.unitize();
					
					Rotation Q = new Rotation(new Vector3D(randomAxis.toArray()), smallAngle);
					for(int i = 0; i < subtreeLeaves.size(); i++){
						int j = subtreeLeaves.get(i);
						Rotation R = new Rotation(new Vector3D(axes[j]), angles[j]);
					
						R = Q.applyTo(R);
						
						axes[j] = R.getAxis().toArray();
						angles[j] = R.getAngle();
					}
					// logProposalRatio is 0 because prior is uniform and proposal is symmetric
					break;
					
				case 1:
					// translation of group of proteins
					System.out.print("Translation: ");
					xlatProposed++;
					double[] shift = new double[3];
					for(int i = 0; i < 3; i++)
						shift[i] = Utils.generator.nextGaussian() * xlatP;
					for(int l = 0; l < subtreeLeaves.size(); l++){
						int j = subtreeLeaves.get(l);
						for(int i = 0; i < 3; i++)
							xlats[j][i] += shift[i];  
					}
					
					// logProposalRatio is 0 because prior is uniform and proposal is symmetric
					
					break;
					
				case 2:
					System.out.print("Library: ");
					Transformation oldSub = new Transformation(axes[index], angles[index], xlats[index]);
					// transformation should be relative to reference protein
					oldSub.xlat = oldSub.xlat.subtract(new ArrayRealVector(xlats[0]));
					Transformation libProp = rotProp.propose(index);
					axes[index] = libProp.rotMatrix.getAxis().toArray();
					angles[index] = libProp.rotMatrix.getAngle();
					xlats[index] = libProp.xlat.toArray();

					// library density 
					logProposalRatio = rotProp.libraryLogDensity(index, oldSub) - 
							rotProp.libraryLogDensity(index, libProp);

					// proposed translation is relative to reference protein
					for(int i = 0; i < 3; i++)
						xlats[index][i] += xlats[0][i];			

					// calculate 'difference' between proposed and current transformations
					double[] diffxlat = new double[3];
					for(int i = 0; i < 3; i++)
						diffxlat[i] = xlats[index][i] - oldxlats[index][i];
					Rotation diffRotMat = libProp.rotMatrix.applyTo(oldSub.rotMatrix.revert());

					for(int i = 0; i < subtreeLeaves.size(); i++){
						int j = subtreeLeaves.get(i);
						if(j != index){
							for(int k = 0; k < 3; k++)
								xlats[j][k] += diffxlat[k];
							Transformation temp = new Transformation(axes[j], angles[j], xlats[j]);
							temp.rotMatrix = diffRotMat.applyTo(temp.rotMatrix);
							axes[j] = temp.rotMatrix.getAxis().toArray();
							angles[j] = temp.rotMatrix.getAngle();
						}
					}
					break;
				}
				
				for(int i = 0; i < subtreeLeaves.size(); i++){
					int j = subtreeLeaves.get(i);
					calcRotation(j);
				}
			}
			
			if (proposingFixedSubtreeAlignment) {
				oldAlign = curAlign;
				logProposalRatio += subtreeRoot.realignToParent();
				curAlign = tree.getState().getLeafAlign();
			}
			
			curLogLike = calcAllColumnContrib();

			if(isParamChangeAccepted(logProposalRatio)) {
				// accepted, nothing to do
				if (param==6) {
					subtreeRotAccept++;
				}
				else if (param==7){
					subtreeAlignAccept++;
				}
				else {
					subtreeRotAlignAccept++;
				}
				System.out.println("accepted!");
			} else {
				// rejected, restore
				curLogLike = oldll;
				int j;
				if (proposingFixedSubtreeRotation) {
					for(int i = 0; i < subtreeLeaves.size(); i++){
						j = subtreeLeaves.get(i);
						axes[j] = oldaxes[j];
						angles[j] = oldangles[j];
						xlats[j] = oldxlats[j];
						rotCoords[j] = oldrots[j];
					}
				}
				if (proposingFixedSubtreeAlignment) {
					subtreeRoot.alignRestore();
					curAlign = oldAlign;
				}
				System.out.println("rejected!");
			}
	
		} else {
			throw new Error("Unknown parameter proposal type");
		}
		
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


