package com.ppfold.algo;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Definition of a node of the phylogenetic tree.
 * 
 * @author Z.Sukosd
 * @see Tree
 */

public class Node implements Serializable {

	private static final long serialVersionUID = 5971946683590961841L;

	// represents a node in the phylogenetic tree
	private double[] Vector;
	private double[] DownBottomVector;
	private double[] DownTopVector;
	private double[] UpBottomVector;
	private double[] UpTopVector;
	
	private double[] tmpvector;
	private double[] tmpvector1;
	private double[][] Matrix;
	private double[] tmpntvector = new double[16];
	private double[][] tmpntmatrix = new double[4][4];

	private double DistanceFromParent;
	private double NewDistanceFromParent; //to store the better distance in MLE calculations
	List<Node> Children = new ArrayList<Node>();
	private String name;
	private int id;

	public Node() {
		Vector = MatrixTools.createVector(4, 1);
	}
	
	public boolean isLeaf() {
		if (Children.size() == 0) {
			return true;
		} else
			return false;
	}

	public void printChildren() {
		for (Node n : this.Children) {
			n.printChildren();
		}
		System.out.format("Id of node: " + id + ",\t name: "
				+ ((name == null || name.equals("")) ? "N/A" : name)
				+ ",\t distance from parent = " + "%.4f" + ",\t has "
				+ Children.size() + " children", DistanceFromParent);
		if (Children.size() > 0) {
			System.out.print("\t with ids: ");
			for (int i = 0; i < this.Children.size(); i++) {
				Node child = Children.get(i);
				System.out.print(child.id
						+ (i == Children.size() - 1 ? "" : ", "));
			}
		}
		System.out.println();

	}

	public int countTotalChildren() {
		int totaloffspringsize = 0;
		for (Node n : this.Children) {
			totaloffspringsize += n.countTotalChildren();
		}
		return this.Children.size() + totaloffspringsize;
	}

	public boolean setChildrenNewBranches(boolean currentvalue) {
		boolean childrensvalue = true;
		for (Node n : this.Children) {
			childrensvalue = childrensvalue&&n.setChildrenNewBranches(currentvalue);
		}
		//currentvalue should be true IF AND ONLY IF:
		//1) this node has finished iterations (ie. value didn't change)
		//2) all its children have finished iterations (ie. n.setChildrenNewBranches() ALWAYS
		//	 returns true)
		currentvalue = (Math.abs(DistanceFromParent-NewDistanceFromParent)<1e-4)&&childrensvalue;
//		if(!currentvalue){
//			System.out.println("Node " + id + " needs more iterations: "
//					+ DistanceFromParent + " vs " + NewDistanceFromParent );
//			System.out.println("Its childrens value " + childrensvalue);
//			System.out.println("Its own value " + (Math.abs(DistanceFromParent-NewDistanceFromParent)<1e-10));
//		}
//		else{
//			System.out.println("Node " + id + " is finished" 
//				+ DistanceFromParent + " vs " + NewDistanceFromParent );
//			System.out.println("Its childrens value " + childrensvalue);
//			System.out.println("Its own value " + (Math.abs(DistanceFromParent-NewDistanceFromParent)<1e-10));
//		
//		}
		this.DistanceFromParent = this.NewDistanceFromParent;
		return currentvalue;
	}
	
	
	
	public void printChildrenMatrix() {
		for (Node n : this.Children) {
			n.printChildrenMatrix();
		}
		System.out.println("Id of node: " + id + ",\t name: "
				+ ((name == null || name.equals("")) ? "N/A" : name)
				+ " the matrix is: ");
		MatrixTools.print(this.Matrix);
	}

	public void calculateChildrenMatrix(double[][] D, double[][] V,
			double[][] V1) {
		for (Node n : this.Children) {
			n.calculateChildrenMatrix(D, V, V1);
		}
		Matrix = calculateMatrix(D, V, V1);

	}

	public double[][] calculateMatrix(double[][] D, double[][] V, double[][] V1) {
		double[][] result = MatrixTools.expRT(D, DistanceFromParent, V, V1);
		return result;
	}

	public void calculateChildrenDownVectors(){
		//this will recursively calculate the vector of all nodes
		//in POSTORDER traversal 
		for(Node n:this.Children){
			n.calculateChildrenDownVectors();
			MatrixTools.resetVector(tmpvector,0);
			MatrixTools.resetVector(tmpvector1,0);
			MatrixTools.copyFromTo(n.DownBottomVector, tmpvector);
			MatrixTools.multiplyMatrixVector(n.Matrix,tmpvector,tmpvector1);
			MatrixTools.multiplySeries(DownBottomVector, tmpvector);
		}
		MatrixTools.copyFromTo(DownBottomVector, DownTopVector);
		MatrixTools.resetVector(tmpvector,0);
		MatrixTools.multiplyMatrixVector(Matrix,DownTopVector,tmpvector);
	}
	
	public void calculateChildrenUpVectors(Node parent) {
		// this will recursively calculate the vector of all nodes
		// in preorder traversal, so first fill in the value, then visit children
		calculateUpVector(parent);
		for (Node n : this.Children) {
			n.calculateChildrenUpVectors(this);
		}
	}

	public void calculateUpVector(Node parent) {
		if(parent!=null){
			//IF NOT ROOT (ie. if it has a parent)
			//up-top is combination of parent's up-bottom....
			MatrixTools.multiplySeries(UpTopVector, parent.UpBottomVector);
			for(Node sibling:parent.getChildren()){
				//...and all siblings' down-top
				if(sibling!=this){
					MatrixTools.multiplySeries(UpTopVector, sibling.DownTopVector);
				}
			}
		}
		//Finally calculate UpBottomVector by doing expRt on UpTop. 
		MatrixTools.copyFromTo(UpTopVector, UpBottomVector);
		MatrixTools.resetVector(tmpvector,0);
		//MatrixTools.multiplyMatrixVector(Matrix,UpBottomVector,tmpvector);
		MatrixTools.multiplyMatrixVector(Matrix,UpBottomVector, tmpvector);
	}
	
	public void calculateChildrenVector() {
		// this will recursively calculate the vector of all nodes
		for (Node n : this.Children) {
			n.calculateChildrenVector();
			// Vector = MatrixTools.multiplySeries(Vector, n.getVector());
			MatrixTools.multiplySeries(Vector, n.getVector());
		}
		calculateVector();
	}

	private void calculateVector() {
		// this will calculate the vector of one node
		MatrixTools.multiplyMatrixVector(Matrix, Vector, tmpvector);
	//	System.out.println("VECTOR of Node " + this.getId() + " is " );
	//	MatrixTools.prints(Vector);
	}

	public void resetChildrenVector(double val) {
		for (Node n : this.Children) {
			n.resetChildrenVector(val);
		}
		for (int i = 0; i < this.Vector.length; i++) {
			Vector[i] = val;
			tmpvector[i] = 0;
		}
	}

	
	public void resetChildrenVector(int m, double val) {
		for (Node n : this.Children) {
			n.resetChildrenVector(m, val);
		}
		double[] vector = new double[m];
		for (int i = 0; i < m; i++) {
			vector[i] = val;
		}
		this.Vector = vector;
		this.tmpvector = new double[m];
	}

	
	public void initializeChildrenUpDownVectors(){
		for (Node n : this.Children) {
			n.initializeChildrenUpDownVectors();
		}
		this.DownBottomVector = new double[4];
		this.DownTopVector = new double[4];
		this.UpBottomVector = new double[4];
		this.UpTopVector = new double[4];
		this.tmpvector = new double[4];
		this.tmpvector1 = new double[4];
	}
	
	public void resetChildrenDownVectors() {
		for (Node n : this.Children) {
			n.resetChildrenDownVectors();
		}
		for (int i = 0; i < 4; i++) {
			DownBottomVector[i] = 1;
			DownTopVector[i] = 1;
			tmpvector[i] = 1;
		}
	}

	public void resetChildrenUpVectors() {
		for (Node n : this.Children) {
			n.resetChildrenUpVectors();
		}
		for (int i = 0; i < 4; i++) {
			UpBottomVector[i] = 1;
			UpTopVector[i] = 1; 
			tmpvector[i] = 1;
		}
	}
	

	public void printChildrenVector() {
		for (Node n : this.Children) {
			n.printChildrenVector();
		}
		System.out.println("Id of node: " + id + ",\t name: "
				+ ((name == null || name.equals("")) ? "N/A" : name)
				+ " the vector is: ");
		MatrixTools.prints(this.Vector);
	}

	public void printChildrenUpDownVectors() {
		for (Node n : this.Children) {
			n.printChildrenUpDownVectors();
		}
		System.out.println("Id of node: " + id + ",\t name: "
				+ ((name == null || name.equals("")) ? "N/A" : name)
				+ " the vectors are (down-bottom, down-top, up-bottom, up-top): ");
		MatrixTools.prints(this.DownBottomVector);
		MatrixTools.prints(this.DownTopVector);
		MatrixTools.prints(this.UpBottomVector);
		MatrixTools.prints(this.UpTopVector);
		MatrixTools.print(this.Matrix);
	}
	
	
	
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Node(double distance) {
		DistanceFromParent = distance;
	}

	public void addChild(Node child) {
		Children.add(child);
	}

	public List<Node> getChildren() {
		return this.Children;
	}

	public double getDistanceFromParent() {
		return DistanceFromParent;
	}

	public void setDistanceFromParent(double distanceFromParent) {
		DistanceFromParent = distanceFromParent;
	}

	public double getNewDistanceFromParent() {
		return NewDistanceFromParent;
	}

	public double[][] getMatrix() {
		return Matrix;
	}

	public void setMatrix(double[][] matrix) {
		Matrix = matrix;
	}

	public double[] getVector() {
		return Vector;
	}

	public void setVector(double[] vector) {
		Vector = vector;
	}

	public double[][] getTmpNtMatrix() {
		return tmpntmatrix;
	}

	public double[] getTmpNtVector() {
		return tmpntvector;
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public void setDownBottomVector(double[] ds) {
		DownBottomVector = ds;		
	}

	public void setUpTopVector(double[] ds) {
		UpTopVector = ds;		
	}
	
	
	public double[] getDownBottomVector(){
		return DownBottomVector;
	}

	public double[] getUpTopVector(){
		return UpTopVector;
	}
	
	public void setDownTopVector(double[][] D, double[][] V, double[][] V1){
		//this sets the DownTopVector to exp(Rt)DownBottomVector, where
		//t = DistanceFromParent. 
		DownTopVector = MatrixTools.multiplyMatrixVector
		(MatrixTools.expRT(D, DistanceFromParent, V, V1), 
				DownBottomVector);
	}
	
	public void setNewDistanceFromParent(double newDistanceFromParent) {
		NewDistanceFromParent = newDistanceFromParent;
	}

	public void setDistanceToOptimum(double distance){
		DistanceFromParent = NewDistanceFromParent;
	}
	
}
