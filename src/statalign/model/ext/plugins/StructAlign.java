package statalign.model.ext.plugins;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

import statalign.base.Tree;
import statalign.base.Vertex;
import statalign.model.ext.ModelExtension;

public class StructAlign extends ModelExtension {

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
	
	/** parameters of structural drift */
	double theta;
	double sigma2;
	
	@Override
	public double logLikeFactor(Tree tree) {
		// simplest (and slowest) approach: get alignment of all leaves and compute likelihood from there
		String[] align = tree.getState().getLeafAlign();
		
		double logli = 0;
		int[] inds = new int[align.length];		// current char indices
		int[] col = new int[align.length];
		for(int i = 0; i < align[0].length(); i++) {
			for(int j = 0; j < align.length; j++)
				col[i] = align[j].charAt(i) == '-' ? -1 : inds[j]++;
			logli += Math.log(columnContrib(col));
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
		for(int i = 0; i < col.length; i++)
			if(col[i] != -1)
				numMatch++;
		// collect indices of ungapped positions
		int[] notgap = new int[numMatch];
		int j = 0;
		for(int i = 0; i < col.length; i++)
			if(col[i] != -1)
				notgap[j++] = i;
		
		// extract covariance corresponding to ungapped positions
		double[][] subCovar = getSubMatrix(fullCovar, notgap, notgap);
		
		// create normal distribution with mean 0 and covariance subCovar
		MultivariateNormalDistribution multiNorm = new MultivariateNormalDistribution(new double[numMatch], subCovar);
	
		double li = 1;
		double[] vals = new double[numMatch];
		// loop over all 3 coordinates
		for(j = 0; j < 3; j++){
			for(int i = 0; i < numMatch; i++)
				vals[i] = rotCoords[notgap[i]][col[notgap[i]]][j];
			li *= multiNorm.density(vals);
		}
		return li;
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
				distMat[i][j] = sigma2 / (2 * theta ) * (1 - Math.exp(2 * theta * distMat[i][j]));
		for(int i = 0; i < tree.names.length; i++)
			for(int j = i + 1; j < tree.names.length; j++)
				distMat[j][i] = distMat[i][j];
		return distMat;
	}
	
	/**
	 * recursive algorithm to traverse tree and calculate distance matrix between leaves 
	 */		
	public int[] calcDistanceMatrix(Vertex vertex, double[][] distMat){
		int[] subTree = new int[distMat.length];
		// either both left and right are null or neither is
		if(vertex.left != null){
			int[] subLeft  = calcDistanceMatrix(vertex.left, distMat);
			int[] subRight = calcDistanceMatrix(vertex.right, distMat);
			int i = 0;
			while(subLeft[i] > -1){
				subTree[i] = subLeft[i];
				i++;
			}
			for(int j = 0; subRight[j] > -1; j++)
				subTree[i+j] = subRight[j];
		}
		else{
			subTree[0] = vertex.index;
			for(int j = 0; j < distMat.length; j++)
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
			for(int j = 0; j < subTree.length; j++)  
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
	public boolean proposeParamChange(Tree tree, double loglike) {
		return false;
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
}
