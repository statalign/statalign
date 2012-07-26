package com.ppfold.algo;

import java.text.DecimalFormat;
import java.util.List;

public class PhyloCalcFuzzy {
	public static double[][] calcSingleColumn(PhyloJobFuzzy phyloJob) {		
		long starttime = System.nanoTime();

		Tree tree = phyloJob.tree;
		List<FuzzyNucleotide[]> columns = phyloJob.columns;
		List<String> names = phyloJob.names;

		double[][] result = new double[columns.size()][1];
		double[][] ntvectors = MatrixTools.createSingleVectors();

		for (int col = 0; col < columns.size(); col++) {
			// Reset all vectors
			tree.getRoot().resetChildrenVector(4, 1);

			// stepping column number.
			// Each column has the nt in the same position from all sequences.
			FuzzyNucleotide [] column = columns.get(col);
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

				// ^a vector of length 4 representing the nucleotide frequencies (frequencies because of ambigious nucleotides)
				double[] vector = node.getVector();
				//MatrixTools.copyFromTo(ntvectors[column[row]], vector);
				// Perhaps should multiple probability by 4, because seems like they give ambiguous nucleotides values of 1 each? Might not matter.
				//double [] frequencies = phyloJob.fuzzyAlignment.getFrequencies(col, row, false);
				//MatrixTools.copyFromTo(frequencies, vector);
				MatrixTools.copyFromTo(column[row].probability, vector);
			}
			// recursively find the vector of all other nodes for this column.
			tree.calculateVectors();

			result[col][0] = MatrixTools.scalarProduct(tree.getRoot()
					.getVector(), phyloJob.param.getPs());
		}
		return result;
	}

	public static double[][] calcDoubleColumn(PhyloJobFuzzy phyloJob) {		
		long starttime = System.nanoTime();
		Tree tree = phyloJob.tree;
		List<FuzzyNucleotide[]> columns = phyloJob.columns;
		List<FuzzyNucleotide[]> columns2 = phyloJob.columns2;
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

				FuzzyNucleotide[] column1 = columns.get(col1);
				FuzzyNucleotide[] column2 = columns2.get(col2);

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
					// ^for nucleotide column1[row]  and nucleotide column2[row]] at position row in column sets create a alphabet*alphabet*16 array, where the vector represents the pairs of nucleotide probabilities, i.e AA, AG, etc.
					// makes a length 16 array, representing paired nucleotide combinations, AA, AC, AG, AT
					double [][] probabilities = new double[column1[row].probability.length][column2[row].probability.length];
					if(phyloJob.fuzzyAlignment.useExpectedFrequencies)
					{
						MatrixTools.multiplyVectorVector(column1[row].probability, column2[row].probability, probabilities);
					}
					else
					{
						probabilities =  phyloJob.fuzzyAlignment.getFrequencyPairs(phyloJob.columnIndices.get(col1), phyloJob.columnIndices2.get(col2), row, true);
					}
										
					//double [][] probabilities =  phyloJob.fuzzyAlignment.getFrequencyPairs(phyloJob.columnIndices.get(col1), phyloJob.columnIndices2.get(col2), row, true);
					
					
					
					double [] serializedProbabilities = new double[16];
					MatrixTools.serializeMatrix(probabilities, serializedProbabilities);
					/*DecimalFormat df = new DecimalFormat("0.000");
					 * System.out.println(col1+ ", "+ col2+":");
					for(int x = 0 ; x < serializedProbabilities.length ; x++)
					{
						System.out.print(df.format(x)+"   ");
					}
					System.out.println();*/
					
					//MatrixTools.copyFromTo(
					//		doublevectors[column1[row]][column2[row]], vector);
					MatrixTools.copyFromTo(serializedProbabilities, vector);
				}
			//System.out.println("$"+phyloJob.columnIndices.get(0)+"#"+phyloJob.columnIndices2.get(0)+"$"+phyloJob.columnIndices.get(phyloJob.columnIndices.size()-1)+"#"+phyloJob.columnIndices2.get(phyloJob.columnIndices2.size()-1));
				//System.out.println();
				//System.out.println(phyloJob.columns.size()+"\t"+phyloJob.columnIndices.size()+"\t"+phyloJob.columns2.size()+"\t"+phyloJob.columnIndices2.size());
				//System.out.println(phyloJob.columnIndices.get(col1)+"\t"+phyloJob.columnIndices2.get(col2));
				//System.out.println(phyloJob.startcol+"\t"+phyloJob.endcol+"\t"+phyloJob.fuzzyAlignment.mapping.size()+"\t"+col1+"\t"+col2+"\t"+(phyloJob.startcol+col1)+"\t"+(phyloJob.startcol+col2));
				// recursively find the vector of all other nodes for this
				// column pair.
				tree.calculateVectors();

				// ^the probability that these columns pair??
				result[col1][col2] = MatrixTools.scalarProduct(
						tree.getRoot().getVector(),
						// ^turns 2 matrix into 1d array
						MatrixTools.serializeMatrix(phyloJob.param.getPd(), tree.getRoot()
								.getTmpNtVector()));
			}
		}
		return result;
	}
}
