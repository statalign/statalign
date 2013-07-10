package statalign.base;

import java.util.Arrays;

import statalign.postprocess.utils.NewickParser;


/**
 * Objects of this class completely describe the state of the MCMC at a certain time.
 * The state space that the chain explores comprises:
 * <ul>
 * <li>a phylogenetic tree (topology + edge lengths)
 * <li>alignments between neighbouring nodes in the tree (implicitly lengths of
 *  the ancestral sequences)
 * <li>parameters of the insertion-deletion model (TKF92) and the user-selected
 *  substitution model
 * </ul>
 * 
 * <p>Apart from the above listed components State objects also store the
 * Felsenstein likelihoods calculated for each character of each sequence and the
 * full likelihood of the state to facilitate postprocessing based on these values
 * (e.g. ancestral sequence estimation, loglikelihood trace etc.).
 * 
 * <p>Tree node representation assumption: first {@link #nl} nodes are the leaves,
 * then come the ancestral nodes.
 * 
 * @author novak
 */
public class State {

	/** Number of nodes */
	public int nn;
	/** Number of leaves (input sequences) */
	public int nl;
	
	/** Root node of the tree */
	public int root;
	
	/** Tree representation: left descendant of each node (first nl are leaves), -1 for none */
	public int[] left;
	/** Tree representation: right descendant of each node (first nl are leaves), -1 for none */
	public int[] right;
	/** Tree representation: parent of each node (first nl are leaves), -1 for none */
	public int[] parent;
	
	/** Edge length to parent, for each node (first nl are leaves) */
	public double[] edgeLen;

	/**
	 * Alignment: for each sequence and each character in the sequence the index of the
	 * ancestral character that it is aligned to, or -i-1 if it is not homologous to
	 * any character in the ancestor (gap) and the next character in the ancestor sequence
	 * (to the right) has index i
	 */
	public int[][] align;
	
	/** Felsenstein likelihoods for each sequence and each character therein */
	public double[][][] felsen;
	/** All sequences as strings, including most likely ancestral characters */
	public String[] seq;
	/** Names of the leaf sequences */
	public String[] name;
	
	/** Indel model parameters (in the order: R, lambda, mu) */
	public double[] indelParams;
	/** Substitution model parameters, as given by the substitution plugin */
	public double[] substParams;
	
	/** Log-likelihood of the state, excluding model extension factors */
	public double logLike;
	
	
	/** Cache to store previously calculated leaf alignment by {@link #getLeafAlign()} */
	private String[] leafAlign;
	/** Cache to store previously calculated full alignment by {@link #getFullAlign()} */
	private String[] fullAlign;
	/** Cache to store previously calculated newick tree string by {@link #getNewickString()} */
	private String newickString;
	/** <code>true</code> if this state comes from the burnin.*/
	public boolean isBurnin;
	
	/**
	 * Constructs a new {@link State} object, filling the given parameter and
	 * pre-allocating arrays (first dimensions only). Does not create parameter arrays.
	 * 
	 * @param nn  number of nodes total (including leaves)
	 */
	public State(int nn) {
		this.nn = nn;
		nl = (nn+1)/2;
		
		left = new int[nn];
		right = new int[nn];
		parent = new int[nn];
		edgeLen = new double[nn];
		
		align = new int[nn][];
		felsen = new double[nn][][];
		seq = new String[nn];
		name = new String[nl];
	}
	
	/**
	 * Returns string representation of alignment between <i>node</i> and its parent. Uses
	 * characters in {@link #seq} to construct the strings, - is printed for gaps.
	 * 
	 * Don't call for root, it will produce array bounds errors.
	 * 
	 * @param node  the node to get the alignment for
	 * @return String  array of two elements, parent comes first
	 */
	public String[] getPairwiseAlign(int node) {
		StringBuilder par = new StringBuilder();	// parent
		StringBuilder des = new StringBuilder();	// descendant
		String s = seq[node], ps = seq[parent[node]];
		
		int[] al = align[node];
		int len = al.length, plen = align[parent[node]].length, i = 0, pi = 0, ali;
		while(i < len || pi < plen) {
			if(i == len || pi < normPos(ali=al[i])) {	// deletion
				par.append(ps.charAt(pi));
				des.append('-');
				pi++;
			} else if(ali < 0) {	// insertion
				par.append('-');
				des.append(s.charAt(i));
				i++;
			} else {	// match
				par.append(ps.charAt(pi));
				des.append(s.charAt(i));
				pi++; i++;
			}
		}
		
		return new String[] { par.toString(), des.toString() };
	}
	
	private int normPos(int pos) {
		return pos < 0 ? -pos-1 : pos;
	}
	
	/**
	 * Returns the multiple alignment of all leaf sequences.
	 */
	public String[] getLeafAlign() {
		if(leafAlign == null) {
			Aligner aligner = new Aligner(true);
			leafAlign = aligner.createAlign();
			
			if(Utils.DEBUG)
				checkConsistency();
		}
		return leafAlign;
	}
	
	/**
	 * Returns the multiple alignment of all sequences, including ancestors.
	 */
	public String[] getFullAlign() {
		if(fullAlign == null) {
			Aligner aligner = new Aligner(false);
			fullAlign = aligner.createAlign();

			if(Utils.DEBUG)
				checkConsistency();
		}
		return fullAlign;
	}
	
	/**
	 * Returns the Newick string representation of the tree
	 */
	public String getNewickString() {
		if(newickString == null) {
			StringBuilder sb = new StringBuilder();
			newick(root, sb);
			newickString = sb.toString();
		}
		return newickString;
	}	
	
	/** Recursively prints newick representation of subtree into sb */
	private void newick(int node, StringBuilder sb) {
		if(node < nl) {		// leaf
            String nameEncoded = NewickParser.getEncodedTaxaName(name[node]);
			sb.append(nameEncoded);
		} else {
			sb.append('(');
			newick(left[node], sb);
			sb.append(',');
			newick(right[node], sb);
			sb.append(')');
		}
		if(node == root) {
			sb.append(';');
		} else {
			sb.append(':');
			sb.append(edgeLen[node]);
		}
	}
	
	private class Aligner {
		
		boolean fullAlign;
		
		char[] column, allGap;
		/** current position in each alignment (indexed by node) */
		int[] pos;
		
		StringBuilder[] rows;
		
		public Aligner(boolean leafOnly) {
			fullAlign = !leafOnly;
			
			int len = leafOnly ? nl : nn;
			allGap = new char[len];
			Arrays.fill(allGap, '-');
			
			rows = new StringBuilder[len];
			for(int i = 0; i < len; i++)
				rows[i] = new StringBuilder();
			
			pos = new int[nn];
		}
		
		String[] createAlign() {
				
			int l = left[root], r = right[root], len = align[root].length, i;
			
			for(i = 0; i < len; i++) {
				before(l, i);
				before(r, i);
				if(fullAlign) {
					newCol();
					column[root] = seq[root].charAt(i);
				}
				at(l, i);
				at(r, i);
				outCol();
			}
			before(l, i);
			before(r, i);
			
			String[] out = new String[rows.length];
			for(i = 0; i < out.length; i++)
				out[i] = rows[i].toString();
			return out;
		}
		
		private void newCol() {
			column = Utils.copyOf(allGap);
		}
		
		private void outCol() {
			if(column == null)
				return;
			for(int i = 0; i < column.length; i++) {
				rows[i].append(column[i]);
			}
			column = null;
		}

		private void before(int node, int pi) {
			int al[] = align[node];
			int i = pos[node], len = al.length, ali = 0, l = left[node], r = right[node];
			boolean leaf = node < nl;
			
			while(i < len && (ali=al[i]) < 0 && -ali-1 == pi) {
				if(!leaf) {
					before(l, i);
					before(r, i);
				}
				if(leaf || fullAlign) {
					newCol();
					column[node] = seq[node].charAt(i);
				}
				if(!leaf) {
					at(l, i);
					at(r, i);
				}
				outCol();
				i++;
			}
			if(!leaf && (
					(i < len && ali >= 0 && ali == pi) ||
					(i == len && pi == align[parent[node]].length))) {
				before(l, i);		// in descendants output insertions (only once!)
				before(r, i);		// before next match and after last column
			}
			pos[node] = i;
		}
		
		private void at(int node, int pi) {
			int al[] = align[node];
			int i = pos[node], len = al.length, ali = 0;
			boolean leaf = node < nl;
			
			if(i < len && (ali=al[i]) >= 0 && ali == pi) {
				if(leaf || fullAlign) {
					if(column == null)
						newCol();
					column[node] = seq[node].charAt(i);
				}
				if(!leaf) {
					at(left[node], i);
					at(right[node], i);
				}
				i++;
			}
			pos[node] = i;
		}
	}
	
	public static void main(String[] args) {
		State s = new State(5);
		s.test();
	}

	private void test() {
		Arrays.fill(parent, -1);
		Arrays.fill(left, -1);
		Arrays.fill(right, -1);
		parent[0] = 3;
		parent[1] = 3;
		parent[2] = 4;
		parent[3] = 4;
		left[3] = 0; right[3] = 1;
		left[4] = 3; right[4] = 2;
		align[0] = new int[] { -2, -2, 1 };
		align[1] = new int[] { 0, 1 };
		align[2] = new int[] { -2, 2 };
		align[3] = new int[] { -1, 2 };
		align[4] = new int[] { 0, 0, 0 };
		root = 4;
		seq[0] = "***";
		seq[1] = "**";
		seq[2] = "**";
		seq[3] = "**";
		seq[4] = "***";
		for(String x : getFullAlign())
			System.out.println(x);
		for(int i = 0; i < 4; i++) {
			System.out.println("///////////");
			for(String x : getPairwiseAlign(i))
				System.out.println(x);
		}
	}

	private void checkConsistency() {
		if(leafAlign == null) {
			Aligner aligner = new Aligner(true);
			leafAlign = aligner.createAlign();
		}
		if(fullAlign == null) {
			Aligner aligner = new Aligner(false);
			fullAlign = aligner.createAlign();
		}
		checkSeqs();
		checkAligns();
	}
	
	private void checkSeqs() {
		// pairwise
		for(int i = 0; i < nn; i++) {
			if(i != root) {
				String[] al = getPairwiseAlign(i);
				checkSeq(al[0], parent[i], "PairP"+parent[i]);
				checkSeq(al[1], i, "PairA"+i);
			}
		}
		
		// leaf align
		for(int i = 0; i < leafAlign.length; i++)
			checkSeq(leafAlign[i], i, "Leaf"+i);
		
		// full align
		for(int i = 0; i < fullAlign.length; i++)
			checkSeq(fullAlign[i], i, "Full"+i);
		
	}

	private void checkSeq(String line, int node, String at) {
		String s = seq[node];
		int p = 0;
		char ch;
		for(int j = 0; j < line.length(); j++) {
			ch = line.charAt(j);
			if(ch == '-')
				continue;
			if(p == s.length() || ch != s.charAt(p))
				throw new Error("Alignment inconsistency id1 at: "+at);
			p++;
		}
		if(p != s.length())
			throw new Error("Alignment inconsistency id2 at: "+at);
	}
	
	private void checkAligns() {
		// full to pairwise
		for(int i = 0; i < nn; i++) {
			if(i != root) {
				String[] pair = getPairwiseAlign(i);
				String[] full = new String[] { fullAlign[parent[i]], fullAlign[i] };
				checkAlign(full, pair, "FullPair"+i);
			}
		}
		// full to leaf
		checkAlign(Arrays.copyOf(fullAlign, leafAlign.length), leafAlign, "FullLeaf");
	}
	
	private void checkAlign(String[] longer, String[] shorter, String at) {
		// remove gaps
		StringBuilder[] sb = new StringBuilder[longer.length];
		int i, j;
		for(i = 0; i < longer.length; i++)
			sb[i] = new StringBuilder();
		for(i = 0; i < longer[0].length(); i++) {
			for(j = 0; j < longer.length; j++)
				if(longer[j].charAt(i) != '-')
					break;
			if(j < longer.length) {
				for(j = 0; j < longer.length; j++)
					sb[j].append(longer[j].charAt(i));
			}
		}
		String[] longNoGaps = new String[longer.length];
		for(i = 0; i < longer.length; i++)
			longNoGaps[i] = sb[i].toString();
		
		// check equality
		for(i = 0; i < longNoGaps.length; i++)
			if(!longNoGaps[i].equals(shorter[i]))
				throw new Error("Alignment inconsistency id3 at: "+at);
	}

	@SuppressWarnings("unused")
	private void print(String[] al) {
		for(String x : al)
			System.out.println(x);
	}
	
}

