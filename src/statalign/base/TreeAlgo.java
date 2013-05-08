package statalign.base;

import java.io.IOException;
import java.io.StringReader;
import java.util.HashMap;

import statalign.base.hmm.HmmNonParam;
import statalign.base.hmm.HmmTkf92;
import statalign.base.thread.Stoppable;
import statalign.base.thread.StoppedException;
import statalign.io.input.plugins.NewickReader;
import statalign.model.subst.SubstitutionModel;

public class TreeAlgo extends Stoppable {
	
	public int GAPOPEN = 9;
	public int GAPEXT = 2;
	public double BIGNUM = 1e100;
	
	public static double MINEDGELEN = 0.01;
	
	public TreeAlgo() {
	}
	
	public TreeAlgo(int gapOpen, int gapExt) {
		this.GAPOPEN = gapOpen;
		this.GAPEXT = gapExt;
	}

	public Tree buildNJTree(String[] sequences, String[] names, SubstitutionModel model, String filename) throws StoppedException {
		Tree tree = new Tree();
		tree.names = names;
		tree.substitutionModel = model;
		tree.title = filename;
		tree.hmm2 = new HmmTkf92(null);
		tree.hmm3 = new HmmNonParam();

		//reading the sequences, transforming them into integer arrays, according to the model
		int[][][] seq = new int[sequences.length][][];
		int[][] which = model.attachedScoringScheme.which;
		for (int i = 0; i < sequences.length; i++) {
			//System.out.println(sequences[i]);
			//System.out.println(ss.which);
			int k = 0;
			for (int j = 0; j < sequences[i].length(); j++) {
				int sum = 0;
				char ch = sequences[i].charAt(j);
				for (int l = 0; l < which[ch].length; l++) {
					sum += which[ch][l];
				}
				if (sum > 0) {
					k++;
				}
			}
			seq[i] = new int[k][model.e.length];
			k = 0;
			for (int j = 0; j < sequences[i].length(); j++) {
				int sum = 0;
				char ch = sequences[i].charAt(j);
				for (int l = 0; l < which[ch].length; l++) {
					sum += which[ch][l];
				}
				if (sum > 0) {
					seq[i][k] = which[ch];
					k++;
				}
			}
			//			for(int j = 0; j < seq[i].length; j++){
			//	System.out.print(seq[i][j]+"\t");
			//}
			//System.out.println();
		}
		// now the pairwise distances
		int[][] dist = new int[seq.length][seq.length];
		int[][] d = null;
		int[][] p = null;
		int[][] q = null;
		int[][] charDist = model.attachedScoringScheme.dist;
		for (int k = 0; k < seq.length; k++) {
			for (int l = 0; l <= k; l++) {
				stoppable();
				if (d == null || d.length < seq[k].length + 1 || d[0].length < seq[l].length + 1) {
					d = new int[seq[k].length + 1][seq[l].length + 1];
					p = new int[seq[k].length + 1][seq[l].length + 1];
					q = new int[seq[k].length + 1][seq[l].length + 1];
				}
				////////////
				d[0][0] = 0;
				for (int i = 1; i <= seq[k].length; i++) {
					d[i][0] = q[i][0] = GAPOPEN - GAPEXT + i * GAPEXT;
				}
				for (int j = 1; j <= seq[l].length; j++) {
					d[0][j] = p[0][j] = GAPOPEN - GAPEXT + j * GAPEXT;
				}
				for (int i = 1; i <= seq[k].length; i++) {
					for (int j = 1; j <= seq[l].length; j++) {
						p[i][j] = Math.min(d[i - 1][j] + GAPOPEN, p[i - 1][j] + GAPEXT);
						q[i][j] = Math.min(d[i][j - 1] + GAPOPEN, q[i][j - 1] + GAPEXT);
						int x, y;
						for (x = 0; seq[k][i - 1][x] == 0; x++) {
							;
						}
						for (y = 0; seq[l][j - 1][y] == 0; y++) {
							;
						}
						d[i][j] = Math.min(d[i - 1][j - 1] + charDist[x][y],
								Math.min(p[i][j], q[i][j]));
					}
				}
				dist[k][l] = dist[l][k] = d[seq[k].length][seq[l].length];
				////////////
				//System.out.println(k+"\t"+l+"\t"+dist[k][l]);
			}

		}

		//// Neighbor Joining algorithm based on the distances calculated above
		// initialization
		tree.vertex = new Vertex[2 * seq.length - 1];
//		for (int i = 0; i < vertex.length; i++) {
//			vertex[i] = new Vertex();
//		}
		double[] sumDist = new double[dist.length];
		for (int i = 0; i < dist.length; i++) {
			sumDist[i] = 0;
			for (int j = 0; j < dist.length; j++) {
				sumDist[i] += dist[i][j];
			}
		}
		int[] where = new int[dist.length];
		for (int i = 0; i < where.length; i++) {
			where[i] = i;
		}
		// the first n vertices will be the leaves
		for (int i = 0; i < seq.length; i++) {
			tree.vertex[i] = new Vertex(tree, i, 0.0, names[i]);
			tree.vertex[i].addSeq(seq[i], sequences[i]);
		}
		// NJ main recursion
		int vnum = seq.length;
		Vertex newVert;
		for (int remN = dist.length; remN > 1; remN--) {
			stoppable();
			double minVal = BIGNUM;
			double val = 0.0;
			int i = -1;
			int j = -1;
			for (int k = 1; k < dist.length; k++) {
				for (int l = 0; l < k; l++) {
					if (where[k] > -1 && where[l] > -1 && (val = (remN - 2) * dist[k][l] - sumDist[k] - sumDist[l]) < minVal) {
						i = k;
						j = l;
						minVal = val;
					}
				}
			}
			newVert = new Vertex(tree, vnum, 0.0, null);	/* new vertex */
			tree.vertex[vnum] = newVert;
			newVert.left = tree.vertex[where[i]];
			newVert.right = tree.vertex[where[j]];
			//System.out.println("Joining vertices "+where[i]+" and "+where[j]);
			newVert.parent = null;
			newVert.left.parent = newVert.right.parent = newVert;
			newVert.left.edgeLength = dist[i][j] / 2 - (remN > 2 ? (sumDist[i] - sumDist[j]) / (2 * remN - 4) : 0.001);
			newVert.right.edgeLength = dist[i][j] - newVert.left.edgeLength;

			val = (newVert.left.getLength() + newVert.right.getLength()) / 0.2;
			newVert.left.edgeLength /= val;
			newVert.right.edgeLength /= val;

			if (newVert.left.edgeLength < MINEDGELEN) {
				newVert.left.edgeLength = MINEDGELEN;
			}
			if (newVert.right.edgeLength < MINEDGELEN) {
				newVert.right.edgeLength = MINEDGELEN;
			}


			//		newVert.left.edgeLength = 0.1;
			//newVert.right.edgeLength = 0.1;

			newVert.left.edgeChangeUpdate();
			newVert.right.edgeChangeUpdate();

			newVert.addEmptySeq();

			newVert.left.fullWin();
			newVert.right.fullWin();
			newVert.fullWin();

			newVert.left.selected = true;
			newVert.right.selected = true;
			newVert.hmm3AlignWithSave();

			//System.out.println("length of the ancestral sequence: "+newVert.length);

			//String[] s = newVert.left.printedAlignment();
			//System.out.println(s[0]+"\n"+s[1]+"\n");

			//s = newVert.right.printedAlignment();
			//System.out.println(s[0]+"\n"+s[1]+"\n");

			where[i] = where[j] = -1;
			sumDist[i] = 0;
			for (int a = 0; a < dist.length; a++) {
				if (where[a] > -1) {
					sumDist[a] -= dist[a][i] + dist[a][j];
					dist[a][i] = dist[i][a] = (dist[a][i] + dist[a][j] - dist[i][j]) / 2;
					sumDist[a] += dist[a][i];
					sumDist[i] += dist[a][i];
				}
			}

			where[i] = vnum;
			vnum++;
		}
		tree.root = tree.vertex[vnum - 1];
		tree.root.calcOrphan();
		tree.root.calcFelsenRecursively();
		
		return tree;
	}
	
	public Tree rearrangeTree(Tree tree, String[] names) {
		// make a copy of the tree
		NewickReader reader = new NewickReader(MINEDGELEN);
		Tree ret;
		try {
			ret = reader.read(new StringReader(tree.printedTree()));
		} catch (IOException ex) {
			return null;
		}

		// rearrange tree's leaf vertices so they come in the input sequence order
		HashMap<String, Integer> nameLookup = new HashMap<String, Integer>();
		int i;
		for(i = 0; i < names.length; i++)
			nameLookup.put(ret.vertex[i].name, i);
		Vertex[] vs = new Vertex[ret.vertex.length];
		for(i = 0; i < names.length; i++) {
			Integer ind = nameLookup.get(names[i]);
			if(ind == null)
				throw new Error("could not find tree node with name "+names[i]);
			vs[i] = ret.vertex[ind];
		}
		for(; i < ret.vertex.length; i++)
			vs[i] = ret.vertex[i];
		ret.vertex = vs;
		ret.names = names;

		return ret;
	}

	/**
	 */
	public void addAlignSeqsToTree(Tree tree, String[] sequences, String[] names, SubstitutionModel model, String filename) throws StoppedException {
		tree.substitutionModel = model;
		tree.title = filename;
		tree.hmm2 = new HmmTkf92(null);
		tree.hmm3 = new HmmNonParam();
		
		// transform sequences into integer arrays, according to the model
		int[][][] seq = new int[sequences.length][][];
		int[][] which = model.attachedScoringScheme.which;
		for (int i = 0; i < sequences.length; i++) {
			//System.out.println(sequences[i]);
			//System.out.println(ss.which);
			int k = 0;
			for (int j = 0; j < sequences[i].length(); j++) {
				int sum = 0;
				char ch = sequences[i].charAt(j);
				for (int l = 0; l < which[ch].length; l++) {
					sum += which[ch][l];
				}
				if (sum > 0) {
					k++;
				}
			}
			seq[i] = new int[k][model.e.length];
			k = 0;
			for (int j = 0; j < sequences[i].length(); j++) {
				int sum = 0;
				char ch = sequences[i].charAt(j);
				for (int l = 0; l < which[ch].length; l++) {
					sum += which[ch][l];
				}
				if (sum > 0) {
					seq[i][k] = which[ch];
					k++;
				}
			}
			//			for(int j = 0; j < seq[i].length; j++){
			//	System.out.print(seq[i][j]+"\t");
			//}
			//System.out.println();
		}
		
		// add leaf sequences to tree
		for(int i = 0; i < names.length; i++) {
			Vertex v = tree.vertex[i];
			v.addSeq(seq[i], sequences[i]);
		}

		// align them recursively
		alignSeqsRec(tree.root);

		tree.root.calcOrphan();
//		tree.root.calcFelsRecursively();
	}

	void alignSeqsRec(Vertex v) throws StoppedException {
		if(v.left == null)
			return;

		// recursively align children first
		alignSeqsRec(v.left);
		alignSeqsRec(v.right);
		
		v.left.edgeChangeUpdate();
		v.right.edgeChangeUpdate();
		
		// initialise ancestor with empty sequence, then do 3-way alignment
		v.addEmptySeq();

		v.left.fullWin();
		v.right.fullWin();
		v.fullWin();

		v.left.selected = true;
		v.right.selected = true;
		stoppable();
		v.hmm3AlignWithSave();
	}

//	int[][] calcDistMatrix(RawSequences seqs, boolean aligned) {
//	}
}
