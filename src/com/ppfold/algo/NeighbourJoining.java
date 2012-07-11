package com.ppfold.algo;

import java.util.ArrayList;
import java.util.List;

/**
 * Contains methods to create a tree by neighbour-joining from the alignment.
 * 
 * @author Z.Sukosd
 */

public class NeighbourJoining {

	static double _PHI = (1 + Math.sqrt(5)) / 2;
	static double _RESPHI = 2 - (1 + Math.sqrt(5)) / 2;
	static int idcounter = -1; //node ID counter
	 
	public static Tree generateTreeNJ(Progress activity, List<String> sequences, List<int[]> columns, 
			List<String> names, Parameters param) throws InterruptedException {
		
		long starttime = System.currentTimeMillis();
		//Number of sequences
		int n = sequences.size();
		System.out.println("Number of sequences: "+ n);
		System.out.println("Generating tree by neighbour-joining...");
		if(n==1){
			//special case, just one sequence
			Node one = new Node(0);
			one.setName(names.get(0));
			return new Tree(one);
		}
		
		//matrix of distances
		double [][] d = new double[n][n]; 
		
		String seq1;
		String seq2;
		activity.setCurrentActivity("Calculating distance matrix...");

		//calculate pairwise distances of sequences
		for(int i=0; i<sequences.size(); i++){
			for(int j=0;j<i+1; j++){
				seq1 = sequences.get(i);
				seq2 = sequences.get(j);
				//System.out.println("Characters: " + seq1 + ", " + seq2);
				activity.checkStop();
				//System.out.println("Iteration " + cnt);
				//create matrices for the branches of the tree.
				double lowerbound = 0;
				double upperbound = 10;
				double midpoint = (upperbound-lowerbound)/2 + _RESPHI * (upperbound - lowerbound);
	
				//System.out.println("Finding distance between " + i + " and " + j);
				d[i][j] = d[j][i] = goldenSectionSearch(activity,
							seq1, seq2, param.getrD(), param.getrV(), 
							param.getrV1(), param.getPr(), 
							lowerbound, midpoint, upperbound, Math.sqrt(1e-4));
				//System.out.println("Distance is found to be " + d[i][j]);
				if(d[i][j] > 10){
					//If the optimum is very large, assume it's infinity and set to 10.
					//(Someone will someday hate me for this comment - sorry ;))
					d[i][j] = d[j][i] = 10;
				}
			}
		}
		
		//System.out.println("Pairwise distance matrix:");
		//MatrixTools.print(d);
		
		//This is to calculate the distances
		//double val = average(d);
		//System.out.println(val + " " + std(d,val) + " " + max(d) + " " + maxdiff(d));
		//Runtime.getRuntime().exit(0);
		
		//do the neighbour-joining  
		activity.setCurrentActivity("Neighbour joining...");
		Tree result = NJ(activity, d,names);
		 
		System.out.println("TOTAL TIME ELAPSED IN NEIGHBOUR-JOINING: " + (System.currentTimeMillis()-starttime)/1000 + " seconds ");
		
//		double value = lnProbability(sequences.get(0), sequences.get(1), 2.3363d, param.getrD(), param.getrV(), 
//				param.getrV1(), param.getPr());
//		System.out.println(value + " = " + Math.exp(value));
//		
//		System.out.println("From A to A: ");
//		double value2 = lnProbability(sequences.get(0), sequences.get(0), 0d, param.getrD(), param.getrV(), 
//				param.getrV1(), param.getPr());
//		System.out.println(value + " = " + Math.exp(value2));
//		
		
		return result;
		
	}

	private static Tree NJ(Progress activity, double [][] d, List<String> names) throws InterruptedException {
		//names contain the names of the leaves. 
		//names.size() corresponds to the number of sequences in the alignment.
		//(each sequence must have a name) 
		int n = names.size();
		
		List<Node> taxa = new ArrayList<Node>();
		//Create a node for each leaf.
		for(int i = 0; i<n; i++){
			Node node = new Node();
			idcounter++;
			node.setId(i);
			//System.out.println(names.get(i) + " has id " + i);
			node.setName(names.get(i));
			node.setDistanceFromParent(0);
			taxa.add(node);
		}
		Node root = joinNeighbourTaxa(activity, taxa,d);
		Tree tree = new Tree(root);
		//tree.print();
		return tree;
	}
	
	private static Node joinNeighbourTaxa(Progress activity, List<Node> taxa, double[][] d) throws InterruptedException{
		activity.checkStop();
		//printTaxa(taxa);
		if(taxa.size()==1){
			//stopping criterion, return the root. 
			if(taxa.get(0).getName()==null){
				//taxa.get(0).setName("root");
			}
			return taxa.remove(0);
		}
		else if(taxa.size()==2){
			//only 2 nodes left.
			//Special case, join them, designate one as root.
			taxa.get(0).addChild(taxa.get(1));
			taxa.get(1).setDistanceFromParent(d[1][0]);
			taxa.get(0).setDistanceFromParent(0);
			taxa.remove(1);
			double[][] newd = new double[1][1];
			newd[0][0] = 1;
			return joinNeighbourTaxa(activity, taxa,newd);		
		}
		else{
			int n = taxa.size(); //how many taxa should be joined
			double partial1 = 0;
			double partial2 = 0;
			double [][] Q = new double [n][n]; 
			for(int i = 0; i<n; i++){
				partial1 = 0;
				for(int k = 0; k<n; k++){
					if(i!=k){
						partial1 = partial1+d[i][k];
					}
				}	
				for(int j = 0; j<i; j++){
					partial2 = 0;
					for(int k = 0; k<n; k++){
						if(j!=k){
						partial2 = partial2+d[j][k];
						}
					}	
			 		Q[i][j] = Q[j][i] = (n-2)*d[i][j] - partial1 - partial2;
			     }
			}
			//System.out.println("Ids:");
			//for(Node tax:taxa){
			//	System.out.print(tax.getId() + " ");
			//}
			//System.out.println();
			//System.out.println("Distance matrix:");
			//MatrixTools.print(d);
			//System.out.println("Q: ");
			//MatrixTools.print(Q);
			
			int[] tojoin = minCoords(Q);
			//System.out.println("Chosen value:" + tojoin[0] + ", " + tojoin[1] + " = " + Q[tojoin[0]][tojoin[1]]);
			//System.out.println("Joining " + taxa.get(tojoin[0]).getId() + "(" + taxa.get(tojoin[0]).getName() + ") " +  
			//		" to " + taxa.get(tojoin[1]).getId()+ "(" + taxa.get(tojoin[1]).getName() + ") ");

			Node newnode = new Node();			
			idcounter++;
			newnode.setId(idcounter);
			//System.out.println("Created node id " + newnode.getId());
			newnode.addChild(taxa.get(tojoin[1]));
			newnode.addChild(taxa.get(tojoin[0]));
			
			//the right coordinates tojoin[1] vs. tojoin[0] are very important here, do not mess! 
			partial1 = 0; //"g" sum, wikipedia definition 
			for(int k = 0; k<n; k++){
				partial1 = partial1+d[tojoin[1]][k];
			}	
			partial2 = 0;//"f" sum, wikipedia definition 
			for(int k = 0; k<n; k++){
				partial2 = partial2+d[tojoin[0]][k];
			}	
			
			//System.out.println("Partial2 (f sum): " + partial2);
			//System.out.println("Partial1 (g sum): " + partial1);
			//System.out.println("Distance of two taxa:" + d[tojoin[0]][tojoin[1]]);
			
			double distance = 0.5*d[tojoin[0]][tojoin[1]] + 
				(partial2-partial1)/(2*n-4); //"f"'s distance from new node
			
			taxa.get(tojoin[0]).setDistanceFromParent(distance);
			//distance of second taxon by reflection
		//	System.out.println(taxa.get(tojoin[0]).getId() + " calculated first ");
			//System.out.println("Distance to new taxon: " + distance);
			double dist1 = distance; 
			
			distance = d[tojoin[0]][tojoin[1]]-distance; //"f"'s distance from the new node
			//System.out.println("Distance of second new taxon: " + distance );
			taxa.get(tojoin[1]).setDistanceFromParent(distance);
			double dist2 = distance; 
			
			//Check for negative branch lengths - in this case just set the negative one to
			//zero and the other one to the whole difference
			if(taxa.get(tojoin[0]).getDistanceFromParent()<0){
				taxa.get(tojoin[0]).setDistanceFromParent(0);
				taxa.get(tojoin[1]).setDistanceFromParent(d[tojoin[0]][tojoin[1]]);
			}
			else if(taxa.get(tojoin[1]).getDistanceFromParent()<0){
				taxa.get(tojoin[1]).setDistanceFromParent(0);
				taxa.get(tojoin[0]).setDistanceFromParent(d[tojoin[0]][tojoin[1]]);	
			}
			else{}
			
			//System.out.println("Children of new node: " + 
			//		newnode.getChildren().get(0).getId() + ":" + newnode.getChildren().get(0).getName()
			//		+ " ( "  + newnode.getChildren().get(0).getDistanceFromParent() + ")" + 
			//		" and " + newnode.getChildren().get(1).getId() + ":" + newnode.getChildren().get(1).getName()
			//		+ " ( "  + newnode.getChildren().get(1).getDistanceFromParent() + ")");
			
			//calculate new distance matrix
			double [][] newd = new double[n-1][n-1];
			//First copy the numbers that are the same. 
			//position in old matrix
			int icnt = 0; 
			int jcnt = 0; 
			for(int i = 0; i<n-2; i++){
				for(int j = 0; j<n-2; j++){
					//icnt=i;
					while(icnt==tojoin[0]||icnt==tojoin[1]){
						//System.out.println("A: Skipping " +icnt);
						icnt++;
					}
					while(jcnt==tojoin[0]||jcnt==tojoin[1]){
						//System.out.println("B: Skipping " + jcnt);
						jcnt++;
					}
					newd[i][j] = d[icnt][jcnt];
				//	System.out.println("Setting " + i + ", " + j + "(" +  newd[i][j] + ") from " + icnt + ", " + jcnt +
				//			"(" +  d[icnt][jcnt] + ")");
				//	System.out.println("NEWD is now: ");
				//	MatrixTools.print(newd);
					jcnt++;
				}
				icnt++;
				jcnt=0;
			}
			//System.out.println("New d before last ones: ");
			//MatrixTools.print(newd);
			
			//Now fill distances to last taxon 
			//(which is always positioned at the end of the new distance table/taxa list) 
			icnt = 0; 
			for(int k = 0; k<n-2; k++){
				//coordinate in new matrix: (i,n-2) and (n-2,i) as it is symmetric.
				while(icnt==tojoin[0]||icnt==tojoin[1]){
					icnt++;
				}
				//System.out.println("Trying to fill " + k + ", " + (n-2) + ", old matrix equivalent " + icnt);
				newd[k][n-2] = newd[n-2][k] = 0.5*(d[tojoin[1]][icnt] + d[tojoin[0]][icnt] - d[tojoin[0]][tojoin[1]]);
				//System.out.println("0.5 * (" + d[tojoin[1]][icnt] + " + " + d[tojoin[0]][icnt] + " - " + d[tojoin[0]][tojoin[1]] + ") = " + newd[k][n-2]);
				icnt++;
			}
			
			//System.out.println("New d:");
			//MatrixTools.print(newd);
			
			//printTaxa(taxa);
			//the second coordinate is always larger so remove that taxon first 
			taxa.remove(tojoin[0]); 
			taxa.remove(tojoin[1]);
			taxa.add(newnode);
			//printTaxa(taxa);
			
			return joinNeighbourTaxa(activity, taxa,newd);
		}
	}
	
	private static void printTaxa(List<Node> taxa){
		System.out.println("TAXON PRINT");
		for(Node taxon:taxa){
			System.out.println("Taxon id: " + taxon.getId() + ", taxon name: " + taxon.getName());
		}
	}
	
	private static int[] minCoords(double [][] matrix){
		//returns the coordinates of the minimum element in the matrix
		if(matrix.length<2){
			System.out.println("Warning: matrix has length 1");
			int [] point = {0, 0};
			return point; 
		}
		double minimum = 100000000; //a very big number. 
		int [] point = {0, 0};  
		for(int i = 0; i<matrix.length; i++){
			for(int j = 0; j<i; j++){
				if(minimum>matrix[i][j]){
					minimum = matrix[i][j];
					point[0] = i; //the second coordinate is always larger 
					point[1] = j;
				}
			}
		}
		return point; 
	}
	
	private static double goldenSectionSearch(Progress act,
			String seq1, String seq2, 
			double[][] D, double[][] V, double [][] V1, double [] Pr,
			double t1, double t2, double t3, double tau) throws InterruptedException{
		// t1 and t3 are the current bounds; the minimum is between them.
		// t2 is the center point, which is closer to t1 than to t3
		
		//System.out.println(lnProbability(seq1, seq2, t2, D, V, V1, Pr));
		//System.out.println("Prob: " + probability(seq1, seq2, t2, D, V, V1, Pr));
		act.checkStop();
		
		// Create a new possible center in the area between t2 and t3, closer to t2
	    double t4 = t2 + _RESPHI * (t3 - t2);
	    
	    if(Math.abs(t3 - t1) < tau * (Math.abs(t2) + Math.abs(t4))){
	        return (t3 + t1) / 2;
	    }

	   // System.out.println(lnProbability(seq1, seq2, t4, D, V, V1, Pr) + 
	   // 		" (for t=" + t4 + ") versus " + lnProbability(seq1, seq2, t2, D, V, V1, Pr)
	   // 		+ " (for t=" + t2 + ")");
	    
	    if(lnProbability(seq1, seq2, t4, D, V, V1, Pr) > 
	       lnProbability(seq1, seq2, t2, D, V, V1, Pr)){
	    //	System.out.println(" A: new bounds "  + t2 + ", " + t3);
	        return goldenSectionSearch(act,seq1, seq2, D, V, V1, Pr, t2, t4, t3, tau);
	    }
	    else{
	    //	System.out.println(" B: new bounds "  + t4 + ", " + t1);
	        return goldenSectionSearch(act,seq1, seq2, D, V, V1, Pr, t4, t2, t1, tau);
	    }
	}
		
	private static double probability(String seq1, String seq2, double t, 
			double[][] D, double[][] V, double [][] V1, double [] Pr){
		
		//Calculate exp{Rt} matrix. 
		double [][] expRT = MatrixTools.expRT(D, t, V, V1);
		double [] vector1 = new double[4];
		double [] vector2 = new double[4];
		double result = 1;
		//Iterate over positions 
		for(int pos = 0; pos<seq1.length(); pos++){
			//Have to take sums of vector positions
			double partial = 0;
			vector1 = MatrixTools.createNtVector(seq1.charAt(pos));
			vector2 = MatrixTools.createNtVector(seq2.charAt(pos));
			for(int vectorpos1 = 0;vectorpos1<4; vectorpos1++){
				if(vector1[vectorpos1]==0){
					continue;
				}
				for(int vectorpos2 = 0;vectorpos2<4;vectorpos2++){
					if(vector2[vectorpos2]==0){
					continue;
				}
				double number = expRT[vectorpos1][vectorpos2]*Pr[vectorpos1]; //partial term;	
				partial = partial+number;
				}
			}
			result = result*partial;
		}
		//System.out.println("Probability: " + result);
		return result; 
	}

	private static double lnProbability(String seq1, String seq2, double t, 
			double[][] D, double[][] V, double [][] V1, double [] Pr){
		
		//Calculate exp{Rt} matrix. 
		double [][] expRT = MatrixTools.expRT(D, t, V, V1);
		double [] vector1 = new double[4];
		double [] vector2 = new double[4];
		double result = 0;
		//Iterate over positions 
		for(int pos = 0; pos<seq1.length(); pos++){
			//Have to take sums of vector positions
			double partial = 0;
			vector1 = MatrixTools.createNtVector(seq1.charAt(pos));
			vector2 = MatrixTools.createNtVector(seq2.charAt(pos));
			for(int vectorpos1 = 0;vectorpos1<4; vectorpos1++){
				if(vector1[vectorpos1]==0){
					continue;
				}
				for(int vectorpos2 = 0;vectorpos2<4;vectorpos2++){
					if(vector2[vectorpos2]==0){
					continue;
				}
				double number = expRT[vectorpos1][vectorpos2]*Pr[vectorpos1]; //partial term;	
				partial = partial+number;
				}
			}
			//System.out.println(partial);
			result = result+Math.log(partial);
		}
		//System.out.println("Log-probability (t = " + t + "): " + result);
		//System.out.println("Probability (t = " + t + "): " + Math.exp(result));
		return result; 
	}
	
	private static double average(double[][] d){
		int n = d.length;
		int cnt = 0; 
		double sum = 0;
		for (int i = 0; i<n; i++){
			for (int j = i+1; j<n; j++){
				sum = sum + d[i][j];
				cnt += 1;
			}
		}
		return sum/cnt;
	}
	
	private static double std(double[][] d, double mean){
		int n = d.length;
		int cnt = 0;
		double sum = 0;
		for (int i = 0; i<n; i++){
			for (int j = i+1; j<n; j++){
				sum = sum + d[i][j]*d[i][j];
				cnt += 1;
			}
		}
		return Math.sqrt(sum/(double)cnt-mean*mean);
	}
	
	private static double max(double[][] d){
		int n = d.length;
		int n2 = n*n;
		double max = 0;
		for (int i = 0; i<n; i++){
			for (int j = i+1; j<n; j++){
				if(d[i][j] > max){
					max = d[i][j];
				}
			}
		}
		return max;
	}
	private static double maxdiff(double[][] d){
		int n = d.length;
		int n2 = n*n;
		double max = 0;
		double min = 10;
		for (int i = 0; i<n; i++){
			for (int j = i+1; j<n; j++){
				if(d[i][j] > max){
					max = d[i][j];
				}
				if(d[i][j] < min){
					min = d[i][j]; 
				}
			}
		}
		return max-min;
	}
	
	
}