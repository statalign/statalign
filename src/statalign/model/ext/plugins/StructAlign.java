package statalign.model.ext.plugins;

import statalign.base.Tree;
import statalign.model.ext.ModelExtension;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

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
	RealMatrix fullCovar;
	
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
		// TODO calc
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
		RealMatrix subCovar = fullCovar.getSubMatrix(notgap, notgap);
		
		// probably need to convert RealMatrix to double[][] first
		multiNorm = MultivariateNormalDistribution(new int[numMatch], subCovar);
		
		double li = 1;
		double[] vals = new double[numMatch];
		// loop over all 3 coordinates
		for(int j = 0; j < 3; j++){
			for(int i = 0; i < numMatch; i++)
				vals[i] = rotCoords[notgap[i]][col[notgap[i]][j];
			li *= multiNorm.density(vals);
		}
		return li;
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
}
