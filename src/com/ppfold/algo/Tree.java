package com.ppfold.algo;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

/**
 * Contains the phylogenetic tree.
 * 
 * @author Z.Sukosd
 * @see Node
 */

public class Tree implements Serializable {

	private static final long serialVersionUID = -5743797481182257647L;
	// the phylogenetic tree
	private Node root;

	private List<Node> leaves;

	public Tree(Node root) {
		this.root = root;
	}

	public void generateLeafList(List<String> names) {
		if(leaves == null){
			//only generate it if it's null 
			leaves = new ArrayList<Node>();
			Node node;
			for (String name : names) {
				node = findSlowlyNodeWithName(name);
				leaves.add(node);
			}
		}
	}

	public Node findNodeWithName(int nr) {
		return leaves.get(nr);
	}

	public Node findSlowlyNodeWithName(String name) {
		Stack<Node> nodelist = new Stack<Node>();
		nodelist.push(root);
		while (!nodelist.isEmpty()) {
			Node node = nodelist.pop();
			if (node.getName() != null && node.getName().equals(name) == true
					&& (
			(node.getName()=="root"&&node.getChildren().size()==0)||
			node.getName()!="root")) {
				//special case for crazy people who want to break my program
				return node;
			} else {
				for (Node n : node.Children) {
					nodelist.push(n);
				}
			}
		}
		return null;
	}

	public List<Node> createListOfNodes(){
		Stack<Node> nodestack = new Stack<Node>();
		List<Node> nodelist = new ArrayList<Node>();
		
		nodestack.push(root);
		while (!nodestack.isEmpty()) {
			Node node = nodestack.pop();
			if (!node.equals(root)) {
				nodelist.add(node);
			}
			for (Node n : node.Children) {
				nodestack.push(n);
			}
		}
		return nodelist;
		
	}
	
	public void calculateVectors() {
		root.calculateChildrenVector();
	}
	
	public void calculateDownVectors(){
		root.calculateChildrenDownVectors();
	}

	public void calculateUpVectors(){
		root.calculateChildrenUpVectors(null);
	}
	
	
	/**
	 * @return the root
	 */
	public Node getRoot() {
		return root;
	}

	/**
	 * @param root
	 *            the root to set
	 */
	public void setRoot(Node root) {
		this.root = root;
	};

	public static Tree copyTree(Tree t) {
		Node root = copyNode(t.getRoot());
		for (Node thisphylonode : t.root.Children) {
			Node child = copyNode(thisphylonode);
			child = copyChildren(thisphylonode, child);
			root.addChild(child);
		}
		Tree newtree = new Tree(root);
		List<String> names = new ArrayList<String>();
		for (Node leaf : t.leaves) {
			names.add(leaf.getName());
		}
		newtree.generateLeafList(names);
		return newtree;
	}

	public static Node copyNode(Node n) {
		// System.out.println("Copying: " + n.getName());
		Node newnode = new Node();
		newnode.setName(n.getName());
		// System.out.println("Name " + n.getName());
		newnode.setId(n.getId());
		newnode.setDistanceFromParent(n.getDistanceFromParent());
		newnode.setMatrix(n.getMatrix());
		return newnode;
	}

	public static Node copyChildren(Node n, Node node) {

		if (n.isLeaf()) {
			return node;
		} else {
			for (Node thisnode : n.Children) {
				Node child = copyNode(thisnode);
				child = copyChildren(thisnode, child);
				node.addChild(child);
			}
			return node;
		}
	}

	public int numberOfNodes() {
		return root.countTotalChildren() + 1;
	}

	public boolean setNewBranches() {
		//argument is "finished iterations"
		return root.setChildrenNewBranches(false);
	}
	
	public void print() {
		root.printChildren();
	}

	public int optimizeBranchLengths(Progress act, List<int[]> columns_int, List<char[]> columns_char, 
			List<String> names, Parameters param, int iterlimit) throws InterruptedException{
		return MaximumLikelihoodTree.optimizeBranchLengths(act, this, columns_int, columns_char, names,param, iterlimit); 
	}
	
}
