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
import org.apache.commons.math3.exception.DimensionMismatchException;
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
import statalign.io.DataType;
import statalign.io.ProteinSkeletons;
import statalign.model.ext.ModelExtManager;
import statalign.model.ext.ModelExtension;
import statalign.postprocess.PluginParameters;

public class StructAlign extends ModelExtension implements ActionListener {
	
	/** the command line identifier of this plugin */
	private static final String CMD_LINE_PLUGIN_ID = "structal";
	
	JToggleButton myButton;

	/** alpha-C atomic coordinate for each sequence and each residue */
	double[][][] coords;
	
	/** alpha-C atomic coordinates under the current set of rotations/translations */
	double[][][] rotCoords;
	
	/** axis of rotation for each sequence */
	double[][] axes;
	/** rotation angle for each protein along the rotation axis */
	double[] rots;
	/** translation vector for each protein */
	double[][] xlats;
	
	/** covariance matrix implied by current tree topology */
	double[][] fullCovar;
	
	/** independence rotation proposal distribution */
	RotationProposal rotProp;
	
	/** parameters of structural drift */
	double theta = .1; // CONSTANT VALUE FOR NOW
	double sigma2 = 1; // SHOULD BE UPDATED WITH MCMC
	
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
	public void init(ModelExtManager manager, PluginParameters params) {
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
	}
	
	@Override
	public double logLikeFactor(Tree tree) {
//		return 0;
		// simplest (and slowest) approach: get alignment of all leaves and compute likelihood from there
		String[] align = tree.getState().getLeafAlign();
		fullCovar = calcFullCovar(tree);
		double logli = 0;
		int[] inds = new int[align.length];		// current char indices
		int[] col = new int[align.length];  // SOMETHING IS WRONG WITH THE WAY INDICES ARE HANDLED HERE
		for(int i = 0; i < align[0].length(); i++) {
			for(int j = 0; j < align.length; j++)
				col[j] = align[j].charAt(i) == '-' ? -1 : inds[j]++;
			logli += columnContrib(col);
		}
		return logli;
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
		if(numMatch == 0)  // CHRIS: this shouldn't happen, but some columns are all gaps
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
		MultiNormCholLog multiNorm = new MultiNormCholLog(new double[numMatch], subCovar);
		
		double logli = 0;
		double[] vals = new double[numMatch];
		// CHRIS: FOR TESTING ONLY, REMOVE WHEN rotCoords PROPERLY FILLED
		rotCoords = coords;
		// loop over all 3 coordinates
		for(j = 0; j < 3; j++){
			for(int i = 0; i < numMatch; i++)
				vals[i] = rotCoords[notgap[i]][col[notgap[i]]][j];
			logli += multiNorm.density(vals);
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
	
	/**
	 * return the full covariance matrix for the tree topology and branch lengths
	 */	
	public double[][] calcFullCovar(Tree tree) {
		// I'm assuming that tree.names.length is equal to the number of vertices here
		double[][] distMat = new double[tree.names.length][tree.names.length];
		calcDistanceMatrix(tree.root, distMat);
		for(int i = 0; i < tree.names.length; i++)
			for(int j = i; j < tree.names.length; j++)
				distMat[i][j] = sigma2 / (2 * theta ) * Math.exp(-theta * distMat[i][j]);
		for(int i = 0; i < tree.names.length; i++)
			for(int j = i + 1; j < tree.names.length; j++)
				distMat[j][i] = distMat[i][j];
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
	
	@Override
	public double proposeParamChange(Tree tree) {
		return 0;
	}

	@Override
	public void afterTreeChange(Tree tree, Vertex nephew, boolean accepted) {
		if(!accepted)
			return;
		// TODO recalculate things
	}
	
	@Override
	public void afterEdgeLenChange(Tree tree, Vertex vertex, boolean accepted) {
		if(!accepted)
			return;
		// TODO recalculate
	}

	/** Adapted from org.apache.commons.math3.distribution.MultivariateNormalDistribution
	 * 
	 * @author Challis
	 * 
	 */
	
	public class MultiNormCholLog{
		/** Dimension. */
		private final int dim;
		/** Vector of means. */
		private final double[] means;
		/** Covariance matrix. */
		private final RealMatrix covarianceMatrix;
		/** The matrix inverse of the covariance matrix. */
		private final RealMatrix covarianceMatrixInverse;
		/** The determinant of the covariance matrix. */
		private final double covarianceMatrixDeterminant;
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
		public MultiNormCholLog(final double[] means,
				final double[][] covariances)
						throws SingularMatrixException,
						DimensionMismatchException,
						NonPositiveDefiniteMatrixException {
			
			this.dim = means.length;

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
			covarianceMatrixDeterminant = covMatDec.getDeterminant();
		}

		/**
		 * Gets the mean vector.
		 *
		 * @return the mean vector.
		 */
		public double[] getMeans() {
			return MathArrays.copyOf(means);
		}

		/**
		 * Gets the covariance matrix.
		 *
		 * @return the covariance matrix.
		 */
		public RealMatrix getCovariances() {
			return covarianceMatrix.copy();
		}

		/** {@inheritDoc} */
		public double density(final double[] vals) throws DimensionMismatchException {
			if (vals.length != dim) {
				throw new DimensionMismatchException(vals.length, dim);
			}

			return -dim / 2 * FastMath.log(2 * FastMath.PI) + 
					-0.5 * FastMath.log(covarianceMatrixDeterminant) +
					getExponentTerm(vals);
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
		
		/* concentration parameters of proposal distribution */
		double kvMF = 500;
		double kvM = 100;
		double s2 = .1;
		
		RotationProposal(){
			// calculate all rotations relative to fixed protein
			libraries = new Transformation[coords.length][];
			for(int i = 0; i < coords.length; i++)
				if(i != fixed)
					libraries[i] = selectBest(calculateAllOptimal(fixed, i, window));
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
					//if(k < 10)
					//	System.out.println(ss);
				}
			}				
			return allOptimal;
		}
		
		public Transformation[] selectBest(Transformation[] library){
			double[] ss = new double[library.length];
			for(int i = 0; i < ss.length; i++)
				ss[i] = library[i].ss;
			double cutoff = new Percentile().evaluate(ss, .01);
			int n = 0;
			for(int i = 0; i < library.length; i++)
				n += library[i].ss < cutoff ? 1 : 0;
			Transformation[] best = new Transformation[n];
			int j = 0;
			for(int i = 0; i < library.length; i++)
				if(library[i].ss < cutoff)
					best[j++] = library[i];
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
		
		/** Calculates the cross product of 2 vectors
		 * 
		 * @param a first vector
		 * @param b second vector
		 * @return cross product
		 */
		RealVector crossprod(RealVector a, RealVector b){
			RealMatrix skew = new Array2DRowRealMatrix(
					new double[][] {{0,a.getEntry(2),-a.getEntry(1)},{-a.getEntry(2),0,a.getEntry(0)},
							{a.getEntry(1),-a.getEntry(0),0}});
			return skew.operate(b);
		}
		
		public class Transformation{
			RealVector axis;
			double rot;
			RealVector xlat;
			RealMatrix rotMatrix;
			// sum of squares for coordinate submatrix
			double ss;
			
			Transformation(){}
			
			Transformation(RealVector axis, double rot){
				this.axis = axis;
				this.rot = rot;
			}
			
			/** Transformation(double[] axis, double rot, double[] xlat){
				this.axis = new ArrayRealVector(axis);
				this.rot = rot;
				this.xlat = new ArrayRealVector(xlat);	
			} */
			

			public Transformation fillAxisAngle(){
				EigenDecomposition eigen = new EigenDecomposition(this.rotMatrix);
				
				double real = 0;
				int i = -1;
				while(real < .99999999){
					i++;
					real = eigen.getRealEigenvalue(i);
				}
				this.axis = eigen.getEigenvector(i);
				
				RealVector v;
				if(this.axis.getEntry(2) < .98)
					v = crossprod(this.axis, new ArrayRealVector(new double[] {0, 0, 1}));
				else
					v = crossprod(this.axis, new ArrayRealVector(new double[] {1, 0, 0}));
				
				RealVector rotv = this.rotMatrix.preMultiply(v);
				
				this.rot = Math.acos(v.dotProduct(rotv));
				this.axis = this.axis.mapMultiply(Math.signum(this.axis.dotProduct(crossprod(v, rotv))));
				return this;
			}
			// Transformation Class
		}
		
		
		/** propose from a library mixture distribution */
		public Transformation propose(Transformation[] library){
			int n = library.length;
			int k = Utils.generator.nextInt(n);
			RealVector axis = simVonMisesFisher(3, kvMF, library[k].axis);
			double rot = simVonMises(kvM, library[k].rot);
			return new Transformation();
		}
		
		public RealVector simVonMisesFisher(int dim, double kappa, RealVector v){
			double b, x0, c, z, u, w;
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
			sample.append(w);
					  
			// calculate axis and angle to rotate (0,0,1) to v
			// axis is cross product of v and (0,0,1)
			RealVector axis = crossprod(v, new ArrayRealVector(new double[]{0,0,1}));
			axis = axis.unitVector();
			
			// shortcut for a.b / |a||b| when a = (0,0,1)
			double angle = Math.acos(v.getEntry(2));
			RealMatrix Q  = new Transformation(axis, angle).fillAxisAngle().rotMatrix;
					  
			// Output simulated vector rotated to center around v
			sample = Q.preMultiply(sample);
			
			return sample;
		}
		
		public double simVonMises(double k2, double theta){
			return 0;
		}
		
		public int accept(){
			// TODO
			return 0;
		}
		
		
		// RotationProposal
	}
	// StructAlign
}


