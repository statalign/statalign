package statalign.postprocess.plugins.contree;

import java.util.ArrayList;
import java.util.BitSet;

import statalign.postprocess.plugins.TreeNode;

public class Cluster extends ArrayList<TreeNode> implements
		Comparable<ArrayList<TreeNode>> {

	// Variables

	public int noOfOccurrences;
	public double edgeLength;

	public BitSet aboveSplit;
	public int nodeRefA;
	public int nodeRefB;
	public boolean added;
	
	// True if is in majority contree...
	public boolean isMajority;

	// Functions

	public Cluster() {
		isMajority = false;
	}

	// Functions

	public int compareTo(ArrayList<TreeNode> other) {
		if (this.size() < other.size())
			return 1;
		if (this.size() > other.size())
			return -1;
		return 0;
	}

}
