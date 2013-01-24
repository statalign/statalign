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

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;

import statalign.base.InputData;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.base.Vertex;
import statalign.base.hmm.Hmm;
import statalign.io.DataType;
import statalign.io.ProteinSkeletons;
import statalign.model.ext.ModelExtension;
import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.PluginParameters;
import statalign.utils.BetaDistribution;
import statalign.utils.GammaDistribution;
import statalign.utils.NormalDistribution;
import cern.jet.math.Bessel;

public class StructAlign extends ModelExtension implements ActionListener {
	
	/** The command line identifier of this plugin */
	private static final String CMD_LINE_PLUGIN_ID = "structal";
	
	JToggleButton myButton;

	/** Alpha-C atomic coordinate for each sequence and each residue */
	double[][][] coords;
	
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
	public boolean globalSigma = false;
	public double epsilon = 5;
	double structTemp = 1;
	
	private int tauInd;
	private int sigma2HInd;
	private int nuInd;
	private int epsilonInd;
	
	
	/** Covariance matrix implied by current tree topology */
	double[][] fullCovar;
	/** Current alignment between all leaf sequences */
	private String[] curAlign;
	/** Current log-likelihood contribution */
	double curLogLike = 0;
	
	private double[][] oldCovar;
	private String[] oldAlign;
	private double oldLogLi;
	
//	public int[] sigProposed;
//	public int[] sigAccept;
//	public int thetaProposed = 0;
//	public int thetaAccept = 0;
	public int sigHProposed;
	public int sigHAccept;
	public int nuProposed;
	public int nuAccept;
	public int tauProposed;
	public int tauAccept;
	public int rotProposed;
	public int rotAccept;
	public int xlatProposed;
	public int xlatAccept;
	public int libProposed;
	public int libAccept;
	public int subtreeRotProposed;
	public int subtreeRotAccept;
	public int subtreeRotAlignProposed;
	public int subtreeRotAlignAccept;
	
	/** independence rotation proposal distribution */
	RotationProposal rotProp;
	
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
	
	/** Constant weights for rotation/translation, sigma2, tau, sigma2Hier, nu, epsilon, subtree rotation and subtree rot+align combined */
	int[] paramPropWConst = { 0, 0, 3, 3, 3, 3, 3, 1 };
	/** Weights per sequence for rotation/translation, sigma2, tau, sigma2Hier, nu, epsilon, subtree rotation  and subtree rot+align combined */
	int[] paramPropWPerSeq = { 5, 3, 0, 0, 0, 0, 0, 0 };
	/** Total weights calculated as const+perseq*nseq */
	int[] paramPropWeights;
	/** Weights for proposing rotation vs translation vs library */
	int[] rotXlatWeights= { 25, 25, 10 };
//	int[] rotXlatWeights= { 25, 25, 0 };	// library off
	
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
	private static final double angleP = 1000;
	// higher values lead to bigger step sizes
	private static final double xlatP = .1;
	
	private static final double MIN_EPSILON = 0.1;
	
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
	public void setActive(boolean active) {
		super.setActive(active);
		System.out.println("StructAlign plugin is now "+(active?"enabled":"disabled"));
	}
	
	@Override
	public void init(PluginParameters params) {
		if(params != null && params.getParameter(CMD_LINE_PLUGIN_ID) != null) {
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
		
		rotProp = new RotationProposal();
		rotCoords = new double[coords.length][][];
		axes = new double[coords.length][];
		angles = new double[coords.length];
		xlats = new double[coords.length][];

		axes[0] = new double[] { 1, 0, 0 };
		angles[0] = 0;
		xlats[0] = new double[] { 0, 0, 0 };
				
		for(i = 1; i < axes.length; i++) {
			Transformation initial = rotProp.propose(i);
			axes[i] = initial.axis.toArray();
			angles[i] = initial.rot;
			xlats[i] = initial.xlat.toArray();
		}
		// alternative initialization
		/*
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
	public double logLikeFactor(Tree tree) {
		String[] align = tree.getState().getLeafAlign();
		checkConsAlign(align);
		curAlign = align;
		
		double[][] covar = calcFullCovar(tree);
		checkConsCovar(covar);
		fullCovar = covar;
		
		if(!checkConsRots())
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
	
	private double calcAllColumnContrib() {
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
		if(numMatch == 0)  // TODO temporary fix because some alignment columns are all gaps
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
	
	private void calcRotation(int ind) {
		double[][] ci = coords[ind], rci = rotCoords[ind];
		if(rci == null)
			rci = rotCoords[ind] = new double[ci.length][];
		Rotation rot = new Rotation(new Vector3D(axes[ind]), angles[ind]);
		for(int i = 0; i < ci.length; i++) {
			rci[i] = rot.applyTo(new Vector3D(ci[i])).add(new Vector3D(xlats[ind])).toArray();
		}
	}

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
		return 25;
	}

	@Override
	public void proposeParamChange(Tree tree) {
		int param = Utils.weightedChoose(paramPropWeights);
		if(param == 0) {
			writeRotationFiles(tree);
			// proposing rotation/translation of a single sequence
			int rotxlat = Utils.weightedChoose(rotXlatWeights);
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
				
				RealMatrix Q = Funcs.calcRotationMatrix(randomAxis, smallAngle);
				RealMatrix R = Funcs.calcRotationMatrix(new ArrayRealVector(axes[ind]), angles[ind]);
				
				R = R.multiply(Q);
				
				Transformation temp = new Transformation();
				temp.rotMatrix = R;
				temp = temp.fillAxisAngle();
				axes[ind] = temp.axis.toArray();
				angles[ind] = temp.rot;
				
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
				axes[ind] = prop.axis.toArray();
				angles[ind] = prop.rot;
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
				if(param == 3)
					sigHAccept++;
					//acceptanceCounts[sigHInd]++;
				else
					nuAccept++;
					//acceptanceCounts[nuInd]++;
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
			
		} else if(param == 6) { // propose rotation to a subtree
			Vertex subtreeRoot = sampleVertex(tree);
			if(subtreeRoot != tree.root){	// if sampleVertex returns the root there are no eligible vertices, skip this step
				subtreeRotProposed++;
				List<Integer> subtreeLeaves = collectLeaves(subtreeRoot);
				int index = subtreeLeaves.get(Utils.generator.nextInt(subtreeLeaves.size()));
			
				double[][] oldaxes = new double[axes.length][axes[0].length];
				double[] oldangles = new double[angles.length];
				double[][] oldxlats = new double[xlats.length][xlats.length];
				double[][][] oldrots = new double[rotCoords.length][rotCoords[0].length][rotCoords[0][0].length];
				int j;
				for(int i = 0; i < subtreeLeaves.size(); i++){
					j = subtreeLeaves.get(i);
					oldaxes[j] = MathArrays.copyOf(axes[j]);
					oldangles[j] = angles[j];
					oldxlats[j] = MathArrays.copyOf(xlats[j]);
					oldrots[j] = rotCoords[j];
					rotCoords[j] = null;	// so that calcRotation creates new array
				}

				double oldll = curLogLike;
				double logProposalRatio = 0;	


				Transformation oldSub = new Transformation(axes[index], angles[index], xlats[index]);
				// transformation should be relative to reference protein
				oldSub.xlat = oldSub.xlat.subtract(new ArrayRealVector(xlats[0]));
				Transformation libProp = rotProp.propose(index);
				axes[index] = libProp.axis.toArray();
				angles[index] = libProp.rot;
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
				oldSub.fillRotationMatrix();
				libProp.fillRotationMatrix();
				RealMatrix diffRotMat = oldSub.rotMatrix.transpose().multiply(libProp.rotMatrix);

				for(int i = 0; i < subtreeLeaves.size(); i++){
					j = subtreeLeaves.get(i);
					if(j != index){
						for(int k = 0; k < 3; k++)
							xlats[j][k] += diffxlat[k];
						Transformation temp = new Transformation(axes[j], angles[j], xlats[j]);
						temp.fillRotationMatrix();
						temp.rotMatrix = temp.rotMatrix.multiply(diffRotMat);
						temp.fillAxisAngle();
						axes[j] = temp.axis.toArray();
						angles[j] = temp.rot;
					}
				}

				for(int i = 0; i < subtreeLeaves.size(); i++){
					j = subtreeLeaves.get(i);
					calcRotation(j);
				}
				curLogLike = calcAllColumnContrib();

				if(isParamChangeAccepted(logProposalRatio)) {
					// accepted, nothing to do
					subtreeRotAccept++;				
				} else {
					// rejected, restore
					for(int i = 0; i < subtreeLeaves.size(); i++){
						j = subtreeLeaves.get(i);
						axes[j] = oldaxes[j];
						angles[j] = oldangles[j];
						xlats[j] = oldxlats[j];
						rotCoords[j] = oldrots[j];
					}
					curLogLike = oldll;
				}
			}
			
		} else if(param == 7) { // propose rotation to a subtree and realignment in a combined step
			System.out.print("Joint proposal: ");
			Vertex subtreeRoot = sampleVertex(tree);
			if(subtreeRoot == tree.root) {		// if sampleVertex returns the root there are no eligible vertices, skip this step
				System.out.println("skipped");
				return;
			}
			subtreeRotAlignProposed++;
			List<Integer> subtreeLeaves = collectLeaves(subtreeRoot);
			int index = subtreeLeaves.get(Utils.generator.nextInt(subtreeLeaves.size()));

			double[][] oldaxes = new double[axes.length][axes[0].length];
			double[] oldangles = new double[angles.length];
			double[][] oldxlats = new double[xlats.length][xlats.length];
			double[][][] oldrots = new double[rotCoords.length][rotCoords[0].length][rotCoords[0][0].length];
			int j;
			for(int i = 0; i < subtreeLeaves.size(); i++){
				j = subtreeLeaves.get(i);
				oldaxes[j] = MathArrays.copyOf(axes[j]);
				oldangles[j] = angles[j];
				oldxlats[j] = MathArrays.copyOf(xlats[j]);
				oldrots[j] = rotCoords[j];
				rotCoords[j] = null;	// so that calcRotation creates new array
			}

			double oldll = curLogLike;
			double logProposalRatio = 0;	


			Transformation oldSub = new Transformation(axes[index], angles[index], xlats[index]);
			// transformation should be relative to reference protein
			oldSub.xlat = oldSub.xlat.subtract(new ArrayRealVector(xlats[0]));
			Transformation libProp = rotProp.propose(index);
			axes[index] = libProp.axis.toArray();
			angles[index] = libProp.rot;
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
			oldSub.fillRotationMatrix();
			libProp.fillRotationMatrix();
			RealMatrix diffRotMat = oldSub.rotMatrix.transpose().multiply(libProp.rotMatrix);

			for(int i = 0; i < subtreeLeaves.size(); i++){
				j = subtreeLeaves.get(i);
				if(j != index){
					for(int k = 0; k < 3; k++)
						xlats[j][k] += diffxlat[k];
					Transformation temp = new Transformation(axes[j], angles[j], xlats[j]);
					temp.fillRotationMatrix();
					temp.rotMatrix = temp.rotMatrix.multiply(diffRotMat);
					temp.fillAxisAngle();
					axes[j] = temp.axis.toArray();
					angles[j] = temp.rot;
				}
			}

			for(int i = 0; i < subtreeLeaves.size(); i++){
				j = subtreeLeaves.get(i);
				calcRotation(j);
			}
//			curLogLike = calcAllColumnContrib();
			
			oldAlign = curAlign;
			
			logProposalRatio += subtreeRoot.realignToParent();
			
			curAlign = tree.getState().getLeafAlign();
			curLogLike = calcAllColumnContrib();

			if(isParamChangeAccepted(logProposalRatio)) {
				// accepted, nothing to do
				subtreeRotAlignAccept++;
				System.out.println("accepted!");
			} else {
				// rejected, restore
				subtreeRoot.alignRestore();
				curAlign = oldAlign;
				
				for(int i = 0; i < subtreeLeaves.size(); i++){
					j = subtreeLeaves.get(i);
					axes[j] = oldaxes[j];
					angles[j] = oldangles[j];
					xlats[j] = oldxlats[j];
					rotCoords[j] = oldrots[j];
				}
				curLogLike = oldll;
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

	/** Adapted from org.apache.commons.math3.distribution.MultivariateNormalDistribution
	 * 
	 * @author Challis
	 * 
	 */
	
	public class MultiNormCholesky{
		/** Dimension. */
		private final int dim;
		/** Vector of means. */
		private final double[] means;
		/** Covariance matrix. */
		private final RealMatrix covarianceMatrix;
		/** The matrix inverse of the covariance matrix. */
		private final RealMatrix covarianceMatrixInverse;
		/** The determinant of the covariance matrix. */
		private double covarianceMatrixDeterminant;
		/** Matrix used in computation of samples. */
		
		/**
		 * Creates a multivariate normal distribution with the given mean vector and
		 * covariance matrix.
		 * <br/>
		 * The number of dimensions is equal to the length of the mean vector
		 * and to the number of rows and columns of the covariance matrix.
		 * It is frequently written as "p" in formulae.
		 *
		 * @param means Vector of means.
		 * @param covariances Covariance matrix.
		 * @throws DimensionMismatchException if the arrays length are
		 * inconsistent.
		 * @throws SingularMatrixException if the eigenvalue decomposition cannot
		 * be performed on the provided covariance matrix.
		 */
		public MultiNormCholesky(final double[] means,
				final double[][] covariances)
						throws SingularMatrixException,
						DimensionMismatchException,
						NonPositiveDefiniteMatrixException {
			
			dim = means.length;

			if (covariances.length != dim) {
				throw new DimensionMismatchException(covariances.length, dim);
			}

			for (int i = 0; i < dim; i++) {
				if (dim != covariances[i].length) {
					throw new DimensionMismatchException(covariances[i].length, dim);
				}
			}

			this.means = MathArrays.copyOf(means);

			covarianceMatrix = new Array2DRowRealMatrix(covariances);

			// Covariance matrix eigen decomposition.
			final CholeskyDecomposition covMatDec;
			try {
				covMatDec = new CholeskyDecomposition(covarianceMatrix);
			}
			catch (NonPositiveDefiniteMatrixException e) {
				System.out.println(e);
				System.out.println("Sigma2 = ");
				for (int i=0; i<sigma2.length; i++) {
					System.out.print(" "+sigma2[i]);
				}
				System.out.println("");
				System.out.println("tau = "+tau);
				System.out.println("epsilon = "+epsilon);
				System.out.println("covariances = ");
				for (int i=0; i<dim; i++) {
					for (int j=0; j<dim; j++) {
						System.out.print(" "+covariances[i][j]);
					}
					System.out.println("");
				}	
			    throw new RuntimeException(e);
			}

			// Compute and store the inverse.
			covarianceMatrixInverse = covMatDec.getSolver().getInverse();
			// Compute and store the determinant.
			covarianceMatrixDeterminant = 0;
			for(int i = 0; i < dim; i++)
				covarianceMatrixDeterminant += 2 * Math.log(covMatDec.getL().getEntry(i, i));			
		}

		/**
		 * Gets the mean vector.
		 *
		 * @return the mean vector.
		 */
		public double[] getMeans() {
			// ADAM: why is this function needed? can't we use means directly instead?
			return MathArrays.copyOf(means);
		}

		/**
		 * Gets the covariance matrix.
		 *
		 * @return the covariance matrix.
		 */
		public RealMatrix getCovariances() {
			// ADAM: same here, but this one is not even used
			return covarianceMatrix.copy();
		}

		/** {@inheritDoc} */
		public double logDensity(final double[] vals) throws DimensionMismatchException {
			if (vals.length != dim) {
				throw new DimensionMismatchException(vals.length, dim);
			}

			double x = (double)-dim / 2 * FastMath.log(2 * FastMath.PI) + 
					// -0.5 * FastMath.log(covarianceMatrixDeterminant) +
					-0.5 * covarianceMatrixDeterminant +
					getExponentTerm(vals);
			
			/*System.out.println(dim);
			System.out.println((double)-dim / 2 * FastMath.log(2 * FastMath.PI));
			System.out.println(-0.5 * covarianceMatrixDeterminant);
			System.out.println(getExponentTerm(vals));*/
			
			return x;
		}

		/**
		 * Computes the term used in the exponent (see definition of the distribution).
		 *
		 * @param values Values at which to compute density.
		 * @return the multiplication factor of density calculations.
		 */
		private double getExponentTerm(final double[] values) {
			final double[] centered = new double[values.length];
			for (int i = 0; i < centered.length; i++) {
				centered[i] = values[i] - getMeans()[i];
			}
			final double[] preMultiplied = covarianceMatrixInverse.preMultiply(centered);
			double sum = 0;
			for (int i = 0; i < preMultiplied.length; i++) {
				sum += preMultiplied[i] * centered[i];
			}
			return -0.5 * sum;
		}
	}
	
	/** Rotation Library
	 * 
	 * @author Challis
	 *
	 */
	
	public class RotationProposal{
		
		/* libraries of axes, angles, and translation vectors for each protein */
		Transformation[][] libraries;
		int window = 25;
		int fixed = 0;
		double percent = 1;  // (value between 0 and 100)
		
		/* concentration parameters of proposal distribution */
		double kvMF = 100;
		double kvM = 100;
		double sd = .1;
		
		// TODO All of the values above this point should probably be chosen elsewhere
		
		RotationProposal(){
			// calculate all rotations relative to fixed protein
			libraries = new Transformation[coords.length][];
			for(int i = 0; i < coords.length; i++)
				if(i != fixed)
					libraries[i] = selectBest(calculateAllOptimal(fixed, i, window), percent);
		}
		
		public Transformation[] calculateAllOptimal(int a, int b, int window){
			RealMatrix A = new Array2DRowRealMatrix(coords[a]);
			RealMatrix B = new Array2DRowRealMatrix(coords[b]);
			int nA = A.getColumn(0).length - window + 1;
			int nB = B.getColumn(0).length - window + 1;
			Transformation[] allOptimal = new Transformation[nA * nB];
			for(int i = 0; i < nA * nB; i++)
				allOptimal[i] = new Transformation();
			
			// create centering matrix 
			RealMatrix C = new Array2DRowRealMatrix(new double[window][window]);
			for(int i = 0; i < window; i++)
				C.setEntry(i, i, 1);
			C = C.scalarAdd(-1/ (double) window);
			
			RealMatrix subA;
			RealMatrix subB;
			
			int k = 0;
			for(int i = 0; i < nA; i++){
				for(int j = 0; j < nB; j++){
					subA = A.getSubMatrix(i, i + window - 1, 0, 2);
					subB = B.getSubMatrix(j, j + window - 1, 0, 2);
					
					SingularValueDecomposition svd = new SingularValueDecomposition(subB.transpose().multiply(C).multiply(subA));
					
					// U V^T from the SVD yields the optimal rotation/reflection matrix
					// for rotations only, use U S V^T, where S is the identity with
					// S_{n,n} = det(U V^T)
					double det = new LUDecomposition(svd.getU().multiply(svd.getVT())).getDeterminant();
					RealMatrix S = new Array2DRowRealMatrix(new double[3][3]);
					S.setEntry(0, 0, 1);
					S.setEntry(1, 1, 1);
					S.setEntry(2, 2, det);
					
					// rotation = U S V^T
					RealMatrix R = svd.getU().multiply(S).multiply(svd.getVT());
					
					// translation = mean(a) - mean(b) R
					RealVector xlat = Funcs.meanVector(subA).subtract(R.preMultiply(Funcs.meanVector(subB)));
					
					// calculate the resulting sum of squares
					double ss = 0;
					
					// a - (R b + 1 xlat)
					RealVector one = new ArrayRealVector(new double[window]).mapAdd(1);
					
					double[][] diff = subA.subtract( subB.multiply(R).add(one.outerProduct(xlat))).getData();
					for(int m = 0; m < diff.length; m++)
						for(int n = 0; n < diff[0].length; n++)
							ss += Math.pow(diff[m][n], 2);
					allOptimal[k].ss = ss;
					allOptimal[k].rotMatrix = R;
					allOptimal[k].xlat = xlat;
					k++;
				}
			}	
			return allOptimal;
		}
		
		/**
		 * @param
		 * @param library - a set of optimal Transformations between protein fragments
		 * @return - the @percent of optimal Transformations with the lowest RMSD
		 */
		
		public Transformation[] selectBest(Transformation[] library, double percent){
			double[] ss = new double[library.length];
			for(int i = 0; i < ss.length; i++)
				ss[i] = library[i].ss;
			double cutoff = new Percentile().evaluate(ss, percent);
			int n = 0;
			for(int i = 0; i < library.length; i++)
				n += library[i].ss < cutoff ? 1 : 0;
			Transformation[] best = new Transformation[n];
			int j = 0;
			for(int i = 0; i < library.length; i++)
				if(library[i].ss < cutoff)
					best[j++] = library[i].fillAxisAngle();
			return best;
		}
		
	
		/** propose from a library mixture distribution */
		public Transformation propose(int index){
			Transformation[] library = libraries[index];
			int n = library.length;
			int k = Utils.generator.nextInt(n);
			Transformation trans = new Transformation();
			trans.axis = vonMisesFisher.simulate(kvMF, library[k].axis);
			trans.rot = vonMises.simulate(kvM, library[k].rot);
			trans.xlat = new ArrayRealVector(3);
			for(int i = 0; i < 3; i++)
				trans.xlat.setEntry(i, Utils.generator.nextGaussian() * sd + library[k].xlat.getEntry(i) );			
			return trans.fillRotationMatrix();
		}
		
		/**
		 * 
		 * @param index - index indicating the set of Transformations on which mixture components are centered
		 * @param candidate - a candidate Transformation whose density is to be evaluated
		 * @return the log density of @candidate according to the library mixture distribution
		 */
		public double libraryLogDensity(int index, Transformation candidate){
			Transformation[] library = libraries[index];
			double density = 0;
			for(int i = 0; i < library.length; i++)
				density += library[i].density(candidate, kvM, kvMF, sd);
			density /= library.length;
			return Math.log(density);
		}						
		// </RotationProposal>
	}
	
	
	static class vonMises{
		/**
		 * 
		 * @param kappa - concentration parameter
		 * @param mean - mean angle
		 * @param angle - angle to be evaluated
		 * @return log density of distribution at @param angle
		 */
		static double logDensity(double kappa, double mean, double angle){
			return kappa * Math.cos(angle - mean) - Math.log(2 * Math.PI * Bessel.i0(kappa));
		}
		
		static double density(double kappa, double mean, double angle){
			return Math.exp(kappa * Math.cos(angle - mean)) / (2 * Math.PI * Bessel.i0(kappa));
		}
		
		/**
		 * 
		 * @param kappa - concentration parameter
		 * @param mean - mean angle
		 * @return - angle simulated from vonMises distribution with given parameters
		 */
		static double simulate(double kappa, double mean){
			double vm = 0, U1, U2, U3, a, b, c, f, r, z;
			boolean cont;
			a = 1 + Math.pow(1 + 4 * Math.pow(kappa, 2), 0.5);
			b = (a - Math.pow(2 * a, 0.5)) / (2 * kappa);
			r = (1 + Math.pow(b, 2)) / (2 * b);
			cont = true;
			while (cont) {
				U1 = Utils.generator.nextDouble();
				z = Math.cos(Math.PI * U1);
				f = (1 + r * z) / (r + z);
				c = kappa * (r - f);
				U2 = Utils.generator.nextDouble();
				if (c * (2 - c) - U2 > 0) {
					U3 = Utils.generator.nextDouble();
					vm = Math.signum(U3 - 0.5) * Math.acos(f) + mean;
					vm = vm % (2 * Math.PI);
					cont = false;
				}
				else {
					if (Math.log(c/U2) + 1 - c >= 0) {
						U3 = Utils.generator.nextDouble();
						vm = Math.signum(U3 - 0.5) * Math.acos(f) + mean;
						vm = vm % (2 * Math.PI);
						cont = false;
					}
				}
			}
			return vm;
		}
		// </vonMises>
	}
	
	static class vonMisesFisher{
		/**
		 * 
		 * @param kappa - concentration parameter
		 * @param mean - mean direction
		 * @param v - vector to be evaluated
		 * @return log density of distribution at @param v
		 */
		static double logDensity(double kappa, RealVector mean, RealVector v){
			return Math.log(kappa) - Math.log(4 * Math.PI * Math.sinh(kappa)) + kappa * mean.dotProduct(v);
		}	
		
		static double density(double kappa, RealVector mean, RealVector v){
			return kappa / (4 * Math.PI * Math.sinh(kappa)) * Math.exp(kappa * mean.dotProduct(v));
		}	
		
		/**
		 * 
		 * @param kappa - concentration parameter
		 * @param v - mean vector
		 * @return - vector simulated from von Mises-Fisher distribution with given parameters
		 */
		// TODO replace with exact simulation method for S2 from Kent et al 2012
		static RealVector simulate(double kappa, RealVector mean){
			double b, x0, c, z, u, w;
			int dim = mean.getDimension();
			RealVector unitSphere = new ArrayRealVector(new double[dim-1]);
			RealVector sample;
			do {
				b = (-2*kappa + Math.pow(Math.pow(4*Math.pow(kappa,2) + (dim-1),2),.5) ) / (dim-1);
				x0 = (1-b)/(1+b);
				c = kappa*x0 + (dim-1)*Math.log(1-Math.pow(x0, 2));
				z = new BetaDistribution((dim-1)/2,(dim-1)/2).sample();
				u = Utils.generator.nextDouble();
				w = (1 - (1+b)*z)/(1 - (1-b)*z);
			} while(kappa*w + (dim-1) * Math.log(1-x0*w) - c < Math.log(u));

			for(int i = 0; i < dim - 1; i++)
				unitSphere.setEntry(i, Utils.generator.nextGaussian());
			unitSphere = unitSphere.unitVector();
			sample = unitSphere.mapMultiply(Math.pow(1-Math.pow(w,2), .5));
			sample = sample.append(w);
					  
			// calculate axis and angle to rotate (0,0,1) to mean
			// axis is cross product of mean and (0,0,1)
			RealVector axis = Funcs.crossProduct(mean, new ArrayRealVector(new double[]{0,0,1}));
			axis = axis.unitVector();
			
			// shortcut for dot(a,b) / |a||b| when a = (0,0,1) and both are unit vectors
			double angle = Math.acos(mean.getEntry(2));
			RealMatrix Q  = Funcs.calcRotationMatrix(axis, angle);
					  
			// Output simulated vector rotated to center around mean
			sample = Q.preMultiply(sample);
			
			return sample;
		}
		// </vonMisesFisher>
	}
	
	/**
	 * 
	 * @author Challis
	 *
	 */
	public class Transformation{
		RealVector axis;
		double rot;
		RealVector xlat;
		RealMatrix rotMatrix;
		// sum of squares for coordinate sub-matrix optimal rotation
		double ss;
		
		Transformation(){}
		
		Transformation(RealVector axis, double rot){
			this.axis = axis;
			this.rot = rot;
		}
		
		Transformation(double[] axis, double rot, double[] xlat){
			this.axis = new ArrayRealVector(axis);
			this.rot = rot;
			this.xlat = new ArrayRealVector(xlat);
		}
		
		/**
		 * 
		 * @return same Transformation object with axis and angle filled
		 */		
		public Transformation fillAxisAngle(){
			EigenDecomposition eigen = new EigenDecomposition(rotMatrix);
			
			double real = 0;
			int i = -1;
			// TODO seems like there should be a better way to do this
			while(real < .99999999){
				i++;
				real = eigen.getRealEigenvalue(i);
			}
			axis = eigen.getEigenvector(i);
			
			RealVector v;
			if(axis.getEntry(2) < .98)
				v = Funcs.crossProduct(axis, new ArrayRealVector(new double[] {0, 0, 1}));
			else
				v = Funcs.crossProduct(axis, new ArrayRealVector(new double[] {1, 0, 0}));
			
			RealVector rotv = rotMatrix.preMultiply(v);
			
			rot = Math.acos(v.dotProduct(rotv) / (v.getNorm() * rotv.getNorm()) );
			// check axis
			axis = axis.mapMultiply(Math.signum(axis.dotProduct(Funcs.crossProduct(v, rotv))));
			return this;
		}
		
		/** 
		 * Calculates the rotation matrix from an axis and angle
		 * with axis x and angle r, the rotation matrix is
		 * cos(r)I + sin(r)cross(x) + (1-cos(r))outer(x)
		 * where cross(x) is the crossproduct matrix of x and outer(x) is the outer (tensor) product
		 * 
		 * @return Transformation object with rotation matrix calculated from axis and angle
		 */
		public Transformation fillRotationMatrix(){
			RealMatrix outer = axis.outerProduct(axis);
			// use transpose of cross product matrix (as given on wikipedia) because
			// we post-mulitply by rotation matrix (other components of final step are symmetric)
			RealMatrix crossTranspose = new Array2DRowRealMatrix(new double[][] { 
					{0, -axis.getEntry(2), axis.getEntry(1)},
					{axis.getEntry(2), 0, -axis.getEntry(0)},
					{-axis.getEntry(1), axis.getEntry(0), 0} 
			});
			RealMatrix Icos = new Array2DRowRealMatrix(new double[3][3]);
			for(int i = 0; i < 3; i++)
				Icos.setEntry(i, i, Math.cos(rot));
			crossTranspose = crossTranspose.scalarMultiply(Math.sin(rot));
			outer = outer.scalarMultiply(1 - Math.cos(rot));
			rotMatrix = Icos.add(crossTranspose).add(outer);
			return this;
		}
		
		public double density(Transformation candidate, double kvM, double kvMF, double sd){
			double density = vonMises.density(kvM, rot, candidate.rot);
			density *= vonMisesFisher.density(kvMF, axis, candidate.axis);
			for(int i = 0; i < 3; i++)
				density *= new NormalDistribution(xlat.getEntry(i), sd).density(candidate.xlat.getEntry(i));
			return density;
		}
		
		// </Transformation>
	}
	
	/**
	 * For storing helpful functions
	 * @author Challis
	 *
	 */
	
	static class Funcs{

		/** Calculates the cross product of 2 vectors
		 * 
		 * @param a first vector
		 * @param b second vector
		 * @return cross product
		 */
		static RealVector crossProduct(RealVector a, RealVector b){
			RealMatrix skew = new Array2DRowRealMatrix(
					new double[][] {{0, -a.getEntry(2), a.getEntry(1)}, {a.getEntry(2), 0, -a.getEntry(0)},
							{-a.getEntry(1), a.getEntry(0), 0}});
			return skew.operate(b);
		}		

		/** 
		 * Calculates the rotation matrix from an axis and angle
		 * with axis x and angle r, the rotation matrix is
		 * cos(r)I + sin(r)cross(x) + (1-cos(r))outer(x)
		 * where cross(x) is the crossproduct matrix of x and outer(x) is the outer (tensor) product
		 * 
		 * @return Transformation object with rotation matrix calculated from axis and angle
		 */
		static RealMatrix calcRotationMatrix(RealVector axis, double rot){
			RealMatrix outer = axis.outerProduct(axis);
			// use transpose of cross product matrix (as given on wikipedia) because
			// we post-multiply by rotation matrix (other components of final step are symmetric)
			RealMatrix crossTranspose = new Array2DRowRealMatrix(new double[][] { 
					{0, axis.getEntry(2), -axis.getEntry(1)},
					{-axis.getEntry(2), 0, axis.getEntry(0)},
					{axis.getEntry(1), -axis.getEntry(0), 0} 
			});
			RealMatrix Icos = new Array2DRowRealMatrix(new double[3][3]);
			for(int i = 0; i < 3; i++)
				Icos.setEntry(i, i, Math.cos(rot));
			crossTranspose = crossTranspose.scalarMultiply(Math.sin(rot));
			outer = outer.scalarMultiply(1 - Math.cos(rot));
			return Icos.add(crossTranspose).add(outer);
		}

		/**
		 * For an n X 3 coordinate matrix, calculate the 1 X 3 mean vector
		 * @param A - coordinate matrix
		 * @return mean vector
		 */
		
		static RealVector meanVector(RealMatrix A){
			RealVector mean = new ArrayRealVector(new double[3]);
			for(int i = 0; i < 3; i ++){
				for(int j = 0; j < A.getColumn(0).length; j++)
					mean.addToEntry(i, A.getEntry(j, i));
				mean.setEntry(i, mean.getEntry(i) / A.getColumn(0).length);
			}
			return mean;
		}
	}

	public Vertex sampleVertex(Tree tree){
		Vertex v;
		int n = tree.vertex.length;		// number of vertices
		int l = coords.length;			// number of leaves
		
		List<Integer> inds = findRefSubtrees(tree, 0);	// returns indices of all ancestor vertices of reference protein
		if(inds.size() + l < n){
			int prop = Utils.generator.nextInt(n-l) + l;		// don't choose a leaf vertex			
			while(inds.contains(new Integer(prop)))
				prop = Utils.generator.nextInt(n-l) + l;		// don't choose a subtree containing reference protein
			v = tree.vertex[prop];
			return v;
		} else
			return tree.root;
	}
	
	public List<Integer> findRefSubtrees(Tree tree, int refInd){
		List<Integer> inds = new ArrayList<Integer>(1);
		moveUp(tree.vertex[refInd].parent, inds);
		return inds;
	}
	
	public void moveUp(Vertex v, List<Integer> inds){
		inds.add(v.index);
		if(v.parent != null)
			moveUp(v.parent, inds);
	}
	
	public List<Integer> collectLeaves(Vertex v){
		List<Integer> inds = new ArrayList<Integer>(1);
		moveDown(v, inds);
		return inds;
	}
	
	public void moveDown(Vertex v, List<Integer> inds){
		if(v.left != null){
			moveDown(v.left, inds);
			moveDown(v.right, inds);
		}
		else
			inds.add(v.index);
	}
	
	public void writeRotationFiles(Tree tree){
		try{
			FileWriter fstream = new FileWriter("names.txt");
			BufferedWriter out = new BufferedWriter(fstream);
			for(int i = 0; i < coords.length; i++)
				out.write(tree.vertex[i].name + "\n");
			out.close();

			FileWriter fstream2 = new FileWriter("trans.txt");
			BufferedWriter out2 = new BufferedWriter(fstream2);
			for(int i = 0; i < coords.length; i++){
				for(int j = 0; j < 3; j++)
					out2.write(xlats[i][j] + "\t");
				out2.write("\n");
			}
			out2.close();
			
			FileWriter fstream3 = new FileWriter("rots.txt");
			BufferedWriter out3 = new BufferedWriter(fstream3);
			for(int i = 0; i < coords.length; i++){
				Transformation trans = new Transformation(axes[i], angles[i], xlats[i]);
				trans.fillRotationMatrix();
				for(int j = 0; j < 3; j++){
					for(int k = 0; k < 3; k++)
						out3.write(trans.rotMatrix.getEntry(j, k) + "\t");
					out3.write("\n");
				}
			}
			out3.close();
			
		} catch (Exception e){
			System.err.println("File writing error: " + e.getMessage());
		}
		
	}
	
	// </StructAlign>
}


