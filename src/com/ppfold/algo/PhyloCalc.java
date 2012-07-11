package com.ppfold.algo;

import java.util.List;

/**
 * Static class, contains calculation instructions of phylogenetic variables.
 * 
 * @author Z.Sukosd
 */

public class PhyloCalc {

	public static double[][] calcSingleColumn(PhyloJob phyloJob) {
		long starttime = System.nanoTime();

		Tree tree = phyloJob.tree;
		List<int[]> columns = phyloJob.columns;
		List<String> names = phyloJob.names;

		double[][] result = new double[columns.size()][1];
		double[][] ntvectors = MatrixTools.createSingleVectors();

		for (int col = 0; col < columns.size(); col++) {
			// Reset all vectors
			tree.getRoot().resetChildrenVector(4, 1);

			// stepping column number.
			// Each column has the nt in the same position from all sequences.
			int[] column = columns.get(col);
			for (int row = 0; row < column.length; row++) {
				// find node corresponding to the rownumber (sequence)
				Node node = tree.findNodeWithName(row);
				if (node == null) {
					System.out.println("Can't find node with name "
							+ names.get(row));
					return null;
				}
				// now set the vector of this Node.
				// the vector of all other nodes will be [1 1 1 1] by default
				// (Node constructor)
				double[] vector = node.getVector();
				MatrixTools.copyFromTo(ntvectors[column[row]], vector);
			}
			// recursively find the vector of all other nodes for this column.
			tree.calculateVectors();

			result[col][0] = MatrixTools.scalarProduct(tree.getRoot()
					.getVector(), phyloJob.param.getPs());
		}
		return result;
	}

	public static double[][] calcDoubleColumn(PhyloJob phyloJob) {
		long starttime = System.nanoTime();
		Tree tree = phyloJob.tree;
		List<int[]> columns = phyloJob.columns;
		List<int[]> columns2 = phyloJob.columns2;
		List<String> names = phyloJob.names;
		int length = columns.size();
		int width = columns2.size();

		double[][] result = new double[length][width];
		double[][][] doublevectors = MatrixTools.createDoubleVectors();
		tree.getRoot().resetChildrenVector(16, 1);
		for (int col1 = 0; col1 < length; col1++) {
			for (int col2 = col1; col2 < width; col2++) {

				// have to reset all vectors.
				tree.getRoot().resetChildrenVector(1);

				int[] column1 = columns.get(col1);
				int[] column2 = columns2.get(col2);

				for (int row = 0; row < column1.length; row++) {
					// find node corresponding to the row number (sequence)
					Node node = tree.findNodeWithName(row);

					if (node == null) {
						System.out.println("Can't find node with name "
								+ names.get(row));
						return null;
					}

					// now set the matrix of this Node.
					// the vector of all other nodes will be [1 1 1 1] by
					// default (Node constructor)
					double[] vector = node.getVector();
					MatrixTools.copyFromTo(
							doublevectors[column1[row]][column2[row]], vector);
				}

				// recursively find the vector of all other nodes for this
				// column pair.
				tree.calculateVectors();

				result[col1][col2] = MatrixTools.scalarProduct(
						tree.getRoot().getVector(), 
						MatrixTools.serializeMatrix(phyloJob.param.getPd(), tree.getRoot()
								.getTmpNtVector()));
			}
		}
		return result;
	}

}
