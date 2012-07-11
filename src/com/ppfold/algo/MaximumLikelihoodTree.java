package com.ppfold.algo;

import java.util.List;


/**
 * Contains methods for the optimization of branch lengths to obtain the maximum 
 * likelihood estimate tree. 
 * 
 * @author Z.Sukosd
 */


public class MaximumLikelihoodTree {
	
	static double _PHI = (1 + Math.sqrt(5)) / 2;
	static double _RESPHI = 2 - (1 + Math.sqrt(5)) / 2;
	
	public static int optimizeBranchLengths(Progress act, Tree tree,List<int[]> columns_int, List<char[]> columns,
			List<String> names,Parameters param, int iterlimit) throws InterruptedException{
		//optimizes the branch lengths of the input tree.
		
		System.out.println("Optimizing branch lengths...");
		long starttime = System.currentTimeMillis();
		double [][] D = param.getrD();
		double [][] V = param.getrV();
		double [][] V1 = param.getrV1();
		double [] Pr = param.getPr();
		
		double [][] upvectors = new double [columns.size()][4];
		double [][] downvectors = new double [columns.size()][4];
		
		boolean iteratemore = true;
		int cnt = 0;

		//System.out.println("Start tree: ");
		//tree.print();
		
		//Create a list of nodes of the tree. 
		List <Node> allnodes = tree.createListOfNodes();
		//Create a list of leaves so they can be referenced quickly later. 
		tree.generateLeafList(names);
		
		double probability = 0;
		
		while(iteratemore&&cnt<iterlimit){
			act.checkStop();
			act.setProgress((double)cnt*(double)1/(double)iterlimit);
			//System.out.println("Iteration " + cnt);
			act.setCurrentActivity("Optimizing tree: iteration " + (cnt+1) + "/" + iterlimit);
			//create matrices for the branches of the tree.
			//these will be the same for each iteration
			//but have to be recalculated after adjusted branch lengths.
			tree.getRoot().calculateChildrenMatrix(D,V,V1);
			tree.getRoot().initializeChildrenUpDownVectors();

			
			for(Node branch_node:allnodes){
				//System.out.println("Optimizing branch length from node " + branch_node.getId() + 
				//		" which is " + branch_node.getDistanceFromParent());

				act.checkStop();
				for (int col = 0; col < columns.size(); col++) {
					// stepping column number.
					// Each column has the nt in the same position from all sequences.
					char[] column = columns.get(col);

					// Reset all vectors for the new column. 
					tree.getRoot().resetChildrenDownVectors();
					tree.getRoot().resetChildrenUpVectors();

					for (int row = 0; row < column.length; row++) {
						// find node corresponding to the rownumber (sequence)
						Node node = tree.findNodeWithName(row);
						if (node == null) {
							System.err.println("Can't find node with name "
									+ names.get(row));
						}
						// now set the down-bottom vectors of this Node.
						node.setDownBottomVector(MatrixTools.createNtVector(column[row]));
						// now set the up-top vector of this node's children
						for(Node n: node.getChildren()){
							n.setUpTopVector(MatrixTools.createNtVector(column[row]));
						}
					}
					// When finished with the whole sequence, 
					// recursively find the down- and up-vectors of all 
					//other nodes for this column.
					tree.calculateDownVectors();
					MatrixTools.copyFromTo(branch_node.getDownBottomVector(), downvectors[col]);
					tree.calculateUpVectors();
					MatrixTools.copyFromTo(branch_node.getUpTopVector(), upvectors[col]);

					//tree.getRoot().printChildrenUpDownVectors();
				}

				//optimize final probability of tree

				probability = lnProbability(upvectors, downvectors, 
						branch_node.getDistanceFromParent(), D, V, V1, Pr);
				//System.out.println("Probability of tree: " + probability);
				
				double lowerbound = 0;
				double upperbound = 10;
				double midpoint = (upperbound-lowerbound)/2 + _RESPHI * (upperbound - lowerbound);
				
				//Change bounds as appropriate 
				branch_node.setNewDistanceFromParent(goldenSectionSearch(
						upvectors, downvectors, param.getrD(), param.getrV(), 
						param.getrV1(), param.getPr(), 
						lowerbound, midpoint, upperbound, Math.sqrt(1e-4)));

				probability = lnProbability(upvectors, downvectors, 
						branch_node.getNewDistanceFromParent(), D, V, V1, Pr);
			//	System.out.println("NEW Probability of tree: " + probability);
			//	System.out.println("OLD length: " + branch_node.getDistanceFromParent());
			//	System.out.println("NEW length: " + branch_node.getNewDistanceFromParent());

			}
			iteratemore = !tree.setNewBranches();
			if(cnt==0){
				System.out.println("Start log-probability of tree: " +  probability);
			}
			cnt++;
		}
		
		if(cnt==iterlimit){
			//System.out.println();
			System.out.println("WARNING! Iteration limit exceeded! (" + iterlimit + ")" + " Tree may not be optimal" +
					" (but it's probably good enough).");
			System.out.println("End log-probability of tree: " + probability);
		}
		else{
			System.out.println("All branch lengths converged after " + cnt + " iterations." );
			System.out.println("End log-probability of tree: " + probability);
		}
		
		//System.out.println(" prob is  = " + Math.exp(probability));
		//System.out.println("Final tree: ");
		//tree.print();
		
		System.out.println("TOTAL TIME ELAPSED IN MLE: " + (System.currentTimeMillis()-starttime)/1000 + " seconds ");

		return cnt;
		
	}
	
	
	private static double goldenSectionSearch(
			double[][] uptopvectors, double[][] downbottomvectors, 
			double[][] D, double[][] V, double [][] V1, double [] Pr,
			double t1, double t2, double t3, double tau){
		// t1 and t3 are the current bounds; the minimum is between them.
		// t2 is the center point, which is closer to t1 than to t3
		
		//System.out.println("Current bounds: " + t1 + ", " + t3 + ", center " + t2);
		
		// Create a new possible center in the area between t2 and t3, closer to t2
	    double t4 = t2 + _RESPHI * (t3 - t2);
	    
	    if(Math.abs(t3 - t1) < tau * (Math.abs(t2) + Math.abs(t4))){
	        return (t3 + t1) / 2;
	    }

	   // System.out.println(probability(uptopvectors, downbottomvectors, t4, D, V, V1, Pr) + 
	   //		" versus " + probability(uptopvectors, downbottomvectors, t2, D, V, V1, Pr));
	    
	    if(lnProbability(uptopvectors, downbottomvectors, t4, D, V, V1, Pr) > 
	       lnProbability(uptopvectors, downbottomvectors, t2, D, V, V1, Pr)){
	    //	System.out.println(" A: new bounds "  + t2 + ", " + t3);
	        return goldenSectionSearch(uptopvectors, downbottomvectors, D, V, V1, Pr, t2, t4, t3, tau);
	    }
	    else{
	    //	System.out.println(" B: new bounds "  + t4 + ", " + t1);
	        return goldenSectionSearch(uptopvectors, downbottomvectors, D, V, V1, Pr, t4, t2, t1, tau);
	    }
	}
		
	/*private static double probability(double[][] uptopvectors, double[][] downbottomvectors, 
		 double t, double[][] D, double[][] V, double [][] V1, double [] Pr){
		//calculates the probability of the tree
		//as the product of  (uptop[i] expRt) .* downbottom[i]  <dot> p_eq, for all columns.

		double probability = 1;
		double number = 0;
		//FOR DEBUG!
		double[] partial = new double[4];
		double[] tmpvector1 = new double[4]; 
		double[][] expRT = MatrixTools.expRT(D, t, V, V1);
		
		//MatrixTools.print(expRT);
		for(int i = 0; i<uptopvectors.length; i++){
			MatrixTools.resetVector(tmpvector1, 0);
			MatrixTools.resetVector(partial, 0);
			MatrixTools.copyFromTo(uptopvectors[i], partial);
			MatrixTools.multiplyVectorMatrix(partial, expRT, tmpvector1);
			MatrixTools.multiplySeries(partial, downbottomvectors[i]);
			number = MatrixTools.scalarProduct(partial, Pr);
			probability = probability*number;
			//System.out.println("Column "+ i + " has probability: " + number);
		}
		//FOR DEBUG:
		return probability;
	}*/
	
	private static double lnProbability(double[][] uptopvectors, double[][] downbottomvectors, 
			 double t, double[][] D, double[][] V, double [][] V1, double [] Pr){
			//calculates the probability of the tree
			//as the product of  (uptop[i] expRt) .* downbottom[i]  <dot> p_eq, for all columns.

			double probability = 0;
			double number = 0;
			//FOR DEBUG!
			double[] partial = new double[4];
			double[] tmpvector1 = new double[4]; 
			double[][] expRT = MatrixTools.expRT(D, t, V, V1);
			//MatrixTools.print(expRT);
			for(int i = 0; i<uptopvectors.length; i++){
				MatrixTools.resetVector(tmpvector1, 0);
				MatrixTools.resetVector(partial, 0);
				MatrixTools.copyFromTo(uptopvectors[i], partial);
				//MatrixTools.copyFromTo(downbottomvectors[i], partial);
				MatrixTools.multiplyVectorMatrix(partial, expRT, tmpvector1);
				MatrixTools.multiplySeries(partial, downbottomvectors[i]);
				//MatrixTools.multiplySeries(partial, uptopvectors[i]);
				number = MatrixTools.scalarProduct(partial, Pr);
				probability = probability+Math.log(number);
				//System.out.println("Column "+ i + " has probability: " + Math.log(number));
			}
			//FOR DEBUG:
			//System.out.println("P(tree | t = " + t + ") = " + probability);
			return probability;
		}
	
	public static int STARTREEoptimizeBranchLengths(Progress act, Tree tree,List<int[]> columns_int, List<char[]> columns,
			List<String> names,Parameters param, int iterlimit) throws InterruptedException{
		//dummy method to set all branch lengths to 0; this is a "hack" to simulate a star-tree
		//for experiments. 
		System.out.println("Setting all branch lengths to zero (simulating star-tree)...");
		List <Node> allnodes = tree.createListOfNodes();
		for(Node branch_node:allnodes){
			branch_node.setDistanceFromParent(0);
		}
		return 0; 
	}
	
}
