package statalign.model.ext.plugins;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JToggleButton;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
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
import cern.jet.math.Bessel;

public class StructAlign extends ModelExtension implements ActionListener {
	
	/** The command line identifier of this plugin */
	private static final String CMD_LINE_PLUGIN_ID = "structal";
	
	JToggleButton myButton;

	/** Alpha-C atomic coordinate for each sequence and each residue */
	double[][][] coords;
	
	/** Alpha-C atomic coordinates under the current set of rotations/translations */
	double[][][] rotCoords;
	
	/** Axis of rotation for each sequence */
	double[][] axes;
	/** Rotation angle for each protein along the rotation axis */
	double[] angles;
	/** Translation vector for each protein */
	double[][] xlats;
	
	/** Covariance matrix implied by current tree topology */
	double[][] fullCovar;
	/** Current alignment between all leaf sequences */
	private String[] curAlign;
	/** Current log-likelihood contribution */
	double curLogLike = 0;
	
	private double[][] oldCovar;
	private String[] oldAlign;
	private double oldLogLi;
	
	int sigProposed = 0;
	int sigAccept = 0;
	int thetaProposed = 0;
	int thetaAccept = 0;
	int rotProposed = 0;
	int rotAccept = 0;
	int xlatProposed = 0;
	int xlatAccept = 0;
	int libProposed = 0;
	int libAccept = 0;
	
	/** independence rotation proposal distribution */
	RotationProposal rotProp;
	
	/** Priors */
	// theta - gamma prior, uses shape/scale parameterization
	GammaDistribution thetaPrior = new GammaDistribution(1, 100);
	
	// sigma2 - gamma prior
	GammaDistribution sigma2Prior = new GammaDistribution(1, 100);
	
	// priors for rotation and translation are uniform
	// so do not need to be included in M-H ratio
	
	/** Proposal tuning parameters */
	// higher values lead to smaller step sizes
	private static final double thetaP = 10;
	private static final double sigma2P = 10;
	private static final double axisP = 100;
	private static final double angleP = 100;
	// higher values lead to bigger step sizes
	private static final double xlatP = .01;
	
	/** Parameters of structural drift */
	public double theta = .1; // TODO CONSTANT VALUE FOR NOW
	public double sigma2 = 1; // SHOULD BE UPDATED WITH MCMC
	
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
			seqMap.put(name, i++);
		coords = new double[inputData.seqs.seqNames.size()][][];
		for(DataType data : inputData.auxData) {
			if(!(data instanceof ProteinSkeletons))
				continue;
			ProteinSkeletons ps = (ProteinSkeletons) data;
			for(i = 0; i < ps.names.size(); i++) {
				String name = ps.names.get(i);
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
		for(i = 0; i < axes.length; i++) {
			axes[i] = new double[] { 1, 0, 0 };
			angles[i] = 0;
			xlats[i] = new double[] { 0, 0, 0 };
		}
	}
	
	@Override
	public double logLikeFactor(Tree tree) {
		String[] align = tree.getState().getLeafAlign();
		checkConsAlign(align);
		curAlign = align;
		
		double[][] covar = calcFullCovar(tree);
		checkConsCover(covar);
		fullCovar = covar;
		
		if(!checkConsRots())
			calcAllRotations();
		
		double logli = calcAllColumnContrib();
		checkConsLogLike(logli);
		curLogLike = logli;
		
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
			logli += columnContrib(col);
		}
		return logli;
	}

	private boolean checkConsAlign(String[] align) {
		if(!Utils.DEBUG || curAlign == null)
			return false;
		if(align.length != curAlign.length)
			throw new Error("Inconsistency in StructAlign, alignment length: "+align.length+", "+curAlign.length);
		for(int i = 0; i < align.length; i++)
			if(!align[i].equals(curAlign[i]))
				throw new Error("Inconsistency in StructAlign, alignment: "+align.length+", "+curAlign.length);
		return true;
	}

	private boolean checkConsCover(double[][] covar) {
		if(!Utils.DEBUG || fullCovar == null)
			return false;
		if(covar.length != fullCovar.length)
			throw new Error("Inconsistency in StructAlign, covar matrix length: "+covar.length+", "+fullCovar.length);
		for(int i = 0; i < covar.length; i++) {
			if(covar[i].length != fullCovar[i].length)
				throw new Error("Inconsistency in StructAlign, covar matrix "+i+" length: "+covar[i].length+", "+fullCovar[i].length);
			for(int j = 0; j < covar[i].length; j++)
				if(Math.abs(covar[i][j]-fullCovar[i][j]) > 1e-5)
					throw new Error("Inconsistency in StructAlign, covar matrix "+i+","+j+" value: "+covar[i][j]+", "+fullCovar[i][j]);
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
		for(j = 0; j < 3; j++){
			for(int i = 0; i < numMatch; i++)
				vals[i] = rotCoords[notgap[i]][col[notgap[i]]][j];
			logli += multiNorm.logDensity(vals);
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
		for(int i = 0; i < tree.names.length; i++)
			for(int j = i; j < tree.names.length; j++)
				distMat[j][i] = distMat[i][j] = sigma2 / (2 * theta ) * Math.exp(-theta * distMat[i][j]);
		return distMat;
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
		addEdgeLength(distMat, subTree, vertex.edgeLength);
		return subTree;
	}
		
	// adds the length of the current edge to the distance between all leaves
	// of a subtree to all other leaves
	// 'rows' contains the indices of vertices in the subtree
	public void addEdgeLength(double[][] distMat, int[] subTree, double edgeLength){
		
		int i = 0;
		while(subTree[i] > -1){
			for(int j = 0; j < distMat.length; j++)  
				distMat[subTree[i]][j] = distMat[subTree[i]][j] + edgeLength;
			i++;		
		}
			
		// edge length should not be added to distance between vertices in the subtree
		// subtract the value from these entries of the distance matrix
		i = 0;
		while(subTree[i] > -1){
			int j = 0;
			while(subTree[j] > -1){
				distMat[subTree[i]][subTree[j]] = 
						distMat[subTree[i]][subTree[j]] - edgeLength;
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
	
	/** Weights for rotation/translation, theta, sigma2, etc. (TODO add all) */
	int[] paramPropWeights = { 10, 3, 3 };
	/** Weights for proposing rotation vs translation vs library */
	int[] rotXlatWeights= { 25, 25, 1 };
//	int[] rotXlatWeights= { 25, 25, 0 };	// library off

	@Override
	public void proposeParamChange(Tree tree) {
		int param = Utils.weightedChoose(paramPropWeights);
		if(param == 0) {
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

			double llratio = 0;

			switch(rotxlat) {
			case 0:
				// rotation of a single sequence
				rotProposed++;
				axes[ind] = vonMisesFisher.simulate(axisP, new ArrayRealVector(axes[ind])).toArray();
				angles[ind] = vonMises.simulate(angleP, angles[ind]);
				
				// llratio is 0 because prior is uniform and proposal is symmetric
				
				break;
			case 1:
				// translation of a single sequence
				xlatProposed++;
				for(int i = 0; i < 3; i++)
					xlats[ind][i] = Utils.generator.nextGaussian() * xlatP + xlats[ind][i];  
				
				// llratio is 0 because prior is uniform and proposal is symmetric
				
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
				llratio = rotProp.libraryLogDensity(ind, old) - 
						  rotProp.libraryLogDensity(ind, prop);
				
				// proposed translation is relative to reference protein
				for(int i = 0; i < 3; i++)
					xlats[ind][i] += xlats[0][i];
							
				break;
			}

			calcRotation(ind);
			curLogLike = calcAllColumnContrib();
			if(isParamChangeAccepted(llratio)) {
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
				if(Utils.DEBUG)
					System.out.println(new String[] { "rot", "xlat", "library" }[rotxlat]+" accepted");
			} else {
				// rejected, restore
				if(Utils.DEBUG)
					System.out.println(new String[] { "rot", "xlat", "library" }[rotxlat]+" rejected");
				axes[ind] = oldax;
				angles[ind] = oldang;
				xlats[ind] = oldxlat;
				rotCoords[ind] = oldrots;
				curLogLike = oldll;
			}
				
		} else {
			
			// proposing new sigma/theta
			double oldpar = param==1 ? theta : sigma2;
			double[][] oldcovar = fullCovar;
			double oldll = curLogLike;
			
			double llratio = 0;
			
			GammaDistribution proposal;
			GammaDistribution reverse;
			
			if(param == 1){
				thetaProposed++;
				// creates a gamma distribution with mean theta & variance controlled by thetaP
				proposal = new GammaDistribution(thetaP, oldpar / thetaP);
				theta = proposal.sample();
				reverse = new GammaDistribution(thetaP, theta / thetaP);
			} else{
				sigProposed++;
				proposal = new GammaDistribution(sigma2P, oldpar / sigma2P);
				sigma2 = proposal.sample();
				reverse = new GammaDistribution(sigma2P, sigma2 / sigma2P);
			}
			
			// TODO do not recalculate distances, only the covariance matrix
			fullCovar = calcFullCovar(tree);
			curLogLike = calcAllColumnContrib();
			
			if(param == 1)
				llratio = Math.log(thetaPrior.density(theta)) + Math.log(reverse.density(oldpar)) 
						  - Math.log(thetaPrior.density(oldpar)) - Math.log(proposal.density(theta));
			else
				llratio = Math.log(sigma2Prior.density(sigma2)) + Math.log(reverse.density(oldpar)) 
				          - Math.log(sigma2Prior.density(oldpar)) - Math.log(proposal.density(sigma2));
			
			if(isParamChangeAccepted(llratio)) {
				if(param == 1)
					thetaAccept++;
				else
					sigAccept++;
				// accepted, nothing to do
				if(Utils.DEBUG)
					System.out.println(new String[] { "theta", "sigma2" }[param-1]+" accepted");
			} else {
				// rejected, restore
				if(Utils.DEBUG)
					System.out.println(new String[] { "theta", "sigma2" }[param-1]+" rejected");
				if(param == 1)
					theta = oldpar;
				else
					sigma2 = oldpar;
				fullCovar = oldcovar;
				curLogLike = oldll;
			}
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
		oldLogLi = curLogLike;
		fullCovar = calcFullCovar(tree);
		curLogLike = calcAllColumnContrib();
		return curLogLike;
	}
	
	@Override
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {
		if(accepted)	// accepted, do nothing
			return;
		// rejected, restore
		fullCovar = oldCovar;
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
			final CholeskyDecomposition covMatDec = new CholeskyDecomposition(covarianceMatrix);

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

			double x = -dim / 2 * FastMath.log(2 * FastMath.PI) + 
//					-0.5 * FastMath.log(covarianceMatrixDeterminant) +
					-0.5 * covarianceMatrixDeterminant +
					getExponentTerm(vals);
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
		int window = 10;
		int fixed = 0;
		double percent = .01;
		
		/* concentration parameters of proposal distribution */
		double kvMF = 500;
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
					RealVector xlat = meanVector(subA).subtract(R.preMultiply(meanVector(subB)));
					
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
		
		/**
		 * For an n X 3 coordinate matrix, calculate the 1 X 3 mean vector
		 * @param A - coordinate matrix
		 * @return mean vector
		 */
		
		public RealVector meanVector(RealMatrix A){
			RealVector mean = new ArrayRealVector(new double[3]);
			for(int i = 0; i < 3; i ++){
				for(int j = 0; j < A.getColumn(0).length; j++)
					mean.addToEntry(i, A.getEntry(j, i));
				mean.setEntry(i, mean.getEntry(i) / A.getColumn(0).length);
			}
			return mean;
		}
	
		/** propose from a library mixture distribution */
		public Transformation propose(int index){
			Transformation[] library = libraries[index];
			int n = library.length;
			int k = Utils.generator.nextInt(n);
			Transformation trans = new Transformation();
			trans.axis = vonMisesFisher.simulate(kvMF, library[k].axis);
			trans.rot = vonMises.simulate(kvM, library[k].rot);
			for(int i = 0; i < 3; i++)
				trans.xlat.setEntry(i, Utils.generator.nextGaussian() * sd + library[k].xlat.getEntry(i));			
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
			
			rot = Math.acos(v.dotProduct(rotv));
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
					{0, axis.getEntry(2), -axis.getEntry(1)},
					{-axis.getEntry(2), 0, axis.getEntry(0)},
					{axis.getEntry(1), -axis.getEntry(0), 0} 
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
					new double[][] {{0,a.getEntry(2),-a.getEntry(1)},{-a.getEntry(2),0,a.getEntry(0)},
							{a.getEntry(1),-a.getEntry(0),0}});
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
	}

	
	// </StructAlign>
}


