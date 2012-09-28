package statalign.base;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import statalign.base.hmm.Hmm2;
import statalign.base.hmm.HmmNonParam;
import statalign.base.hmm.HmmSilent;
import statalign.base.hmm.HmmTkf92;
import statalign.base.thread.Stoppable;
import statalign.base.thread.StoppedException;
import statalign.model.score.SubstitutionScore;
import statalign.model.score.plugins.Blosum62;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Dayhoff;
import statalign.postprocess.plugins.contree.CNetwork;

/**
 * This is the current tree in the MCMC run.
 * It extends Stoppable, so the MCMC run can be terminated in graphical interface mode
 * while it calls a function in the Tree class.
 * @author miklos, novak
 */
public class Tree extends Stoppable {

    /** Characters on this tree undergo substitutions according to this model. */
    public SubstitutionModel substitutionModel;
    /**
     * When two fragments of two sequences on two neighbor vertices of the tree
     * are realigned, the proposed new alignment is drawn from a forward-backward
     * sampling of this HMM. It is a simplified version of the TKF92 model.
     */
    public Hmm2 hmm2;
    /**
     * When two fragments of two sequences on two sibling vertices are realigned
     * drawing a novel ancestral substring for their parent vertex, the proposed new
     * alignment is drawn from a forward-backward sampling of this HMM. It is
     * a pair-HMM having states emitting into the non-observable ancestral sequence.
     */
    public HmmSilent hmm3;
    /** The array of vertices of the tree. */
    public Vertex vertex[];

    /** The root of the tree */
    public Vertex root;

    /** The name of the tree */
    public String title;

    /** The name of the sequences */
    public String[] names;

    static final int GAPOPEN = 9;
    static final int GAPEXT = 2;
    static final double BIGNUM = 1e100;
    
    /** The heat parameter for this MCMC chain. */
	public double heat = 1.0d;

	/* TODO: what the ...? */
    public CNetwork network;


    Tree() {
        try {
            substitutionModel = new Dayhoff();
        } catch (IOException e) {
        }

    }

    /**
     * Constructor to read a tree and an alignment from an output file.
     * Not used in the current version. It aims at letting the possibility to
     * get an input tree and an alignment as the starting state of the MCMC.
     */
    // TODO: before using check if leaf vertices are created first!
    Tree(File f, SubstitutionScore ss) {
        try {
            substitutionModel = new Dayhoff();
            hmm2 = new HmmTkf92(null);
            hmm3 = new HmmNonParam();
            BufferedReader bf = new BufferedReader(new FileReader(f));
            String s = bf.readLine().replaceFirst("\\w+\t", "");
            String[] s1 = s.split("\t");
            hmm2.params[0] = Double.parseDouble(s1[3]);
            hmm2.params[1] = Double.parseDouble(s1[5]);
            hmm2.params[2] = Double.parseDouble(s1[7]);
            s = bf.readLine().replaceFirst("\\w+\t", "");
            Tree t = new Tree(s);
            bf.readLine();                    // skip `Alignment' header
            root = t.root;
            vertex = t.vertex;
            for (int i = 0; i < vertex.length; i++) {
                vertex[i].owner = this;
            }
            for (int i = 0; i < vertex.length; i++) {
                vertex[i].edgeChangeUpdate();
            }
            // reading the alignment
            String[] alignment = new String[vertex.length];
            for (int i = 0; i < alignment.length; i++) {
                alignment[i] = bf.readLine().replaceFirst("\\w+\t", "");
            }
            Vertex actual = root;
            actual.selected = false;
            //    System.out.println("alignment length: "+alignment.length+" Tree: "+printedTree());
            for (int i = 0; i < alignment.length; i++) {
                String temp = alignment[i].substring(31, alignment[i].length());
                //	System.out.println(i+" "+temp);
                AlignColumn column = new AlignColumn(actual);
                column.seq = new double[substitutionModel.e.length];
                actual.first = column;
                actual.length = 0;
                for (int j = 0; j < temp.length(); j++) {
                    if (temp.charAt(j) != '-') {
                        //	System.out.println("character in position: "+temp.charAt(j));
                        //column.seq[ss.which[temp.charAt(j)]] = 1.0;
                        for (int k = 0; k < column.seq.length; k++) {
                            column.seq[k] = ss.which[temp.charAt(j)][k];
                        }
                        AlignColumn next = new AlignColumn(actual);
                        next.seq = new double[substitutionModel.e.length];

                        column.next = next;
                        next.prev = column;
                        column = next;
                        actual.length++;
                    }
                }
                /*
            AlignColumn next = new AlignColumn(actual);
            next.seq = new double[substitutionModel.e.length];
            column.next = next;
            next.prev = column;
            column = next;
            actual.length++;
                     */
                actual.last = column;

                // here we use this integer variable for remembering to the index in the alignment...
                // pretty messy :)
                actual.leafCount = i;
                //alignment on the tree
                if (actual.parent != null) {
                    AlignColumn parent = actual.parent.first;
                    column = actual.first;
                    String temp1 = alignment[actual.parent.leafCount].substring(31);
                    // System.out.println(temp1+"\n"+temp);
                    for (int j = 0; j < temp.length(); j++) {
                        if (temp.charAt(j) != '-' && temp1.charAt(j) != '-') {
                            column.parent = parent;
                            column.orphan = false;
                            if (actual.parent.left == actual) {
                                parent.left = column;
                            } else {
                                parent.right = column;
                            }
                            column = column.next;
                            parent = parent.next;
                        } else if (temp.charAt(j) == '-' && temp1.charAt(j) != '-') {
                            parent = parent.next;
                        } else if (temp.charAt(j) != '-' && temp1.charAt(j) == '-') {
                            column.parent = parent;
                            column.orphan = true;
                            column = column.next;
                        }
                    }
                    column.parent = parent;
                    column.orphan = false;
                    if (actual.parent.left == actual) {
                        parent.left = column;
                    } else {
                        parent.right = column;
                    }

                } else {//we are at root, set all aligncolumns to orphan
                    column = actual.first;
                    while (column != null) {
                        column.orphan = true;
                        column = column.next;
                    }
                }

                actual.selected = false;
                if (actual.left == null && actual.right == null) {
                    actual = actual.parent;
                    while (actual.selected && actual.parent != null) {
                        actual = actual.parent;
                    }
                    actual.selected = true;
                    actual = actual.right;
                } else {
                    actual = actual.left;
                }
            }

            //double check:
            //String[] check = printedAlignment();
            //for(int i = 0; i < check.length; i++){
            //	System.out.println(check[i]);
            //}
            root.calcFelsRecursively();
            root.calcIndelLikeRecursively();
            //System.out.println("Log-likelihood: "+getLogLike());

            names = new String[(vertex.length + 1) / 2];
            int index = 0;
            for (int i = 0; i < vertex.length; i++) {
                if (vertex[i].left == null) {
                    names[index] = vertex[i].name;
                    index++;
                }
            }

        } catch (IOException e) {
        }
    }

    void compareLike(Tree tree) {
        System.out.println("this.logli=" + getLogLike());
        System.out.println("tree.logli=" + tree.getLogLike());
        System.out.println("this.indel=" + root.indelLogLike);
        System.out.println("tree.indel=" + tree.root.indelLogLike);
        System.out.println("this.orphan=" + root.orphanLogLike);
        System.out.println("tree.orphan=" + tree.root.orphanLogLike);
        root.compareLike(tree.root);
    }

    /**
     * This constructor constructs a tree using a newick representation of the tree. It <b>does not</b>
     * construct its HMMs and substitution model!!!
     */
    Tree(String descriptor) {
        int numberofnodes = 1;
        for (int i = 0; i < descriptor.length(); i++) {
            if (descriptor.charAt(i) == ':') {
                numberofnodes++;
            }
        }
        vertex = new Vertex[numberofnodes];
        root = new Vertex();
        root.old = new Vertex();

        String leftDescriptor = "";
        int counter = 0;
        int j;
        for (j = 1; counter != 0 || descriptor.charAt(j) != ','; j++) {
            leftDescriptor += descriptor.charAt(j);
            if (descriptor.charAt(j) == '(') {
                counter++;
            }
            if (descriptor.charAt(j) == ')') {
                counter--;
            }
        }
        //System.out.println(leftDescriptor);
        Vertex leftChild = new Vertex(leftDescriptor, root);
        root.left = leftChild;

        String rightDescriptor = "";
        for (j += 1; counter != -1; j++) {
            rightDescriptor += descriptor.charAt(j);
            if (descriptor.charAt(j) == '(') {
                counter++;
            }
            if (descriptor.charAt(j) == ')') {
                counter--;
            }
        }
        rightDescriptor = rightDescriptor.substring(0, rightDescriptor.length() - 1);
        //System.out.println(rightDescriptor);
        Vertex rightChild = new Vertex(rightDescriptor, root);
        root.right = rightChild;


        //System.out.println(printedTree());
        //Mapping vertices
        int actualNumber = numberofnodes - 1;
        Vertex actualVertex = root;
        root.selected = false;
        vertex[numberofnodes - 1] = root;
        while (actualNumber != 0) {
            //    System.out.println(" we are in the root of: "+actualVertex.print()+" it is selected: "+actualVertex.selected+
            //		       " actualNumber: "+actualNumber);
            if (actualVertex.left == null && actualVertex.right == null) {
                actualVertex.selected = true;
                while (actualVertex.selected) {
                    actualVertex = actualVertex.parent;
                }
                actualVertex.selected = true;
            } else {
                actualVertex = (actualVertex.selected ? actualVertex.right : actualVertex.left);
                actualNumber--;
                vertex[actualNumber] = actualVertex;
                actualVertex.selected = false;
            }
        }

        names = new String[(vertex.length + 1) / 2];
        int index = 0;
        for (int i = 0; i < vertex.length; i++) {
            if (vertex[i].left == null) {
                names[index] = vertex[i].name;
                index++;
            }
        }

    }

    /**
     * This constructor generates a tree and puts aligned sequences onto it.
     * The sequences are aligned using an iterative alignment scheme and Neighbor Joining.
     * First pairwise alignments are used to calculate distances, and these distances are used
     * to construct a tree using Neighbor Joining. As the tree is being constructed, sequences
     * are aligned together in an iterative manner.
     * @param sequences The sequences themselves.
     * @param names     The name of the sequences. The code calling this constructor is responsible to
     *                  have the same order of sequences and their names, as it is assumed that the name of
     *                  sequences[i] is names[i].
     * @param model     The time-continuous Markov model describing the substitution process of
     *                  sequences
     * @param ss        The substitution score matrix that is used to obtain pairwise distances.
     * @param filename  The name of the file that contained the sequences. This will be the name of
     *                  the tree, appearing in the graphical interface showing multiple alignments.
     */
    public Tree(String[] sequences, String[] names, SubstitutionModel model, SubstitutionScore ss, String filename) throws StoppedException {
        this.names = names;
        substitutionModel = model;
        title = filename;
        hmm2 = new HmmTkf92(null);
        hmm3 = new HmmNonParam();

        //reading the sequences, transforming them into integer arrays, according to the model
        int[][][] seq = new int[sequences.length][][];
        for (int i = 0; i < sequences.length; i++) {
            //System.out.println(sequences[i]);
            //System.out.println(ss.which);
            int k = 0;
            for (int j = 0; j < sequences[i].length(); j++) {
                int sum = 0;
                char ch = sequences[i].charAt(j);
                for (int l = 0; l < ss.which[ch].length; l++) {
                    sum += ss.which[ch][l];
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
                for (int l = 0; l < ss.which[ch].length; l++) {
                    sum += ss.which[ch][l];
                }
                if (sum > 0) {
                    seq[i][k] = ss.which[ch];
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
        int[][] charDist = ss.dist;
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
        vertex = new Vertex[2 * seq.length - 1];
//        for (int i = 0; i < vertex.length; i++) {
//            vertex[i] = new Vertex();
//        }
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
            vertex[i] = new Vertex(this, 0.0, seq[i], names[i], sequences[i]);
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
            newVert = new Vertex(this, 0.0);    /* new vertex */
            vertex[vnum] = newVert;
            newVert.left = vertex[where[i]];
            newVert.right = vertex[where[j]];
            //System.out.println("Joining vertices "+where[i]+" and "+where[j]);
            newVert.parent = null;
            newVert.left.parent = newVert.right.parent = newVert;
            newVert.left.edgeLength = dist[i][j] / 2 - (remN > 2 ? (sumDist[i] - sumDist[j]) / (2 * remN - 4) : 0.001);
            newVert.right.edgeLength = dist[i][j] - newVert.left.edgeLength;

            val = (newVert.left.length + newVert.right.length) / 0.2;
            newVert.left.edgeLength /= val;
            newVert.right.edgeLength /= val;

            if (newVert.left.edgeLength < 0.01) {
                newVert.left.edgeLength = 0.01;
            }
            if (newVert.right.edgeLength < 0.01) {
                newVert.right.edgeLength = 0.01;
            }


            //	    newVert.left.edgeLength = 0.1;
            //newVert.right.edgeLength = 0.1;

            newVert.left.edgeChangeUpdate();
            newVert.right.edgeChangeUpdate();

            AlignColumn fake = new AlignColumn(newVert);
            newVert.first = fake;
            newVert.last = fake;
            newVert.length = 0;
            fake.left = newVert.left.last;
            fake.right = newVert.right.last;
            newVert.left.last.parent = fake;
            newVert.left.last.orphan = false;
            newVert.right.last.parent = fake;
            newVert.right.last.orphan = false;

            //			/* FAKE!!!! alignment */
            //			AlignColumn prev = new AlignColumn(newVert);
            //			newVert.first = prev;
            //			prev.seq = new double[model.e.length];
            //			AlignColumn left = newVert.left.first;
            //			AlignColumn right = newVert.right.first;
            //			left.parent = prev;
            //			left.orphan = false;
            //			right.parent = prev;
            //			right.orphan = false;
            //			left = left.next;
            //			right = right.next;
            //			while(left.next != null && right.next != null){
            //				AlignColumn actual = new AlignColumn(newVert);
            //				actual.seq = new double[model.e.length];
            //				actual.prev = prev;
            //				prev.next = actual;
            //				left.parent = actual;
            //				left.orphan = false;
            //				right.parent = actual;
            //				right.orphan = false;
            //				right = right.next;
            //				left = left.next;
            //				prev = actual;
            //			}
            //			AlignColumn actual = new AlignColumn(newVert);
            //			actual.seq = new double[model.e.length];
            //			actual.prev = prev;
            //			prev.next = actual;
            //
            //			newVert.last = actual;

            //??????

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
        root = vertex[vnum - 1];
        root.calcOrphan();
        root.calcFelsRecursively();
        /////////////

    }

    public double getLogLike() {
        return root.indelLogLike + root.orphanLogLike;
    }

    public int countLeaves() {
        return root.countLeaves();
    }


    int getTopVertexId(int ind) {
        if (root == vertex[ind]) {
            return 0;
        }
        if (root.left == vertex[ind]) {
            return 1;
        }
        if (root.right == vertex[ind]) {
            return 2;
        }
        return -1;
    }
    
    /**
     * Returns a {@link State} object representing the current state. Assumes that
     * leaves come first in the {@link #vertex} array.
     */
	public State getState() {
		int nn = vertex.length;
		int nl = (nn+1)/2;
		State state = new State(nn);
		
		int i;
		HashMap<Vertex, Integer> lookup = new HashMap<Vertex, Integer>();
		for(i = 0; i < nn; i++)
			lookup.put(vertex[i], i);
		
		Vertex v;
		state.root = lookup.get(root);
		for(i = 0; i < nn; i++) {
			v = vertex[i];
			state.left[i] = v.left != null ? lookup.get(v.left) : -1;
			state.right[i] = v.right != null ? lookup.get(v.right) : -1;
			state.parent[i] = v.parent != null ? lookup.get(v.parent) : -1;
			state.edgeLen[i] = v.edgeLength;

			state.align[i] = v.getAlign();
			state.felsen[i] = v.getFelsen();
			state.seq[i] = v.sequence();
		}
		for(i = 0; i < nl; i++)
			state.name[i] = vertex[i].name;
		
		state.indelParams = hmm2.params.clone();
		state.substParams = substitutionModel.params.clone();
		state.logLike = getLogLike();
		
		return state;
	}

    /**
     * Generates a String contaning the description of the tree in Newick format.
     * @return the string contaning the description of the tree in Newick format.
     */
    public String printedTree() {
        return root.print();
    }
    
	public double getLogPrior() {
		return - root.calcSumOfEdges() - Math.log(substitutionModel.getPrior()) - 
				hmm2.params[1] - hmm2.params[2];
	}

    /**
     * Returns with a String array containing the alignment of sequences on the tree.
     * @param type This can be
     *             <ul>
     *             <li><tt>StatAlign</tt> Our own format. Contains the ancestral sequences, too
     *             <li><tt>Clustal</tt> Standard Clustal alignment format.
     *             <li><tt>Fasta</tt> Standard fasta format, alignments are broken into 60 character long lines.
     *             <li><tt>Phylip</tt> Standard Phylip format
     *             <li><tt>Nexus</tt> Nexus format. MrBayes can read it, but we did not test it for other programs.
     *             </ul>
     * @return The string array containing the multiple alignment in the prescribed format.
     */
//    public String[] printedAlignment(String type) {
//        String[] s = printedAlignment();
//        return Utils.alignmentTransformation(s, type, this);
//
//    }

    private String[] printedAlignment() {
        String[] s = root.printedMultipleAlignment();
        /*
          int p = 29;
          while(allSpace(s,p)){
              for(int i = 0; i < s.length; i++){
                  String temp = s[i].substring(0, p) + s[i].substring(p+1,s[i].length());
                  s[i] = temp;
              }
              p--;
          }
          p++;
          for(int i = 0; i < s.length; i++){
              s[i] = s[i].substring(0,p)+" "+s[i].substring(p, s[i].length());
          }
          */
        return s;

    }

    /**
     * It generates the array of row sequences.
     * The format is sequence name TAB row sequence.
     * @return The string array containing the sequence names and sequences.
     */
    public String rowSequences() {
        String s = "";
        for (int i = 0; i < vertex.length; i++) {
            if (vertex[i].left == null) {
                s += vertex[i].name + "\t" + vertex[i].sequence() + "\n";
            }
        }
        return s;
    }

    /** This function is only for testing purposes. */
    public static void main(String[] args) {
        /*
      Tree tree = new Tree(new String[] {"kkkkkk", "kkkkkkkk", "kkkkkkkk", "kkkkkkkk"},
                   new String[] {"A","B","C", "D"},
                   new NonParametric(),
                   new Blosum62());
      System.out.println(tree.printedTree());
      tree.root.calcFelsRecursively();
      System.out.println(tree.root.orphanLogLike);
           */
        try {
            Tree tree = new Tree(new File("newick"), new Blosum62());
            System.out.println(tree.printedTree());
        } catch (FileNotFoundException e) {
        } catch (IOException e) {
        }
    }
    
    /** 
     * It roots the tree on a taxon 
     * @param taxonName Taxon on which to root tree 
     */ 
    public void makeRoot(String taxonName){ 
        if(vertex.length>0){ 
                Vertex taxon = vertex[0]; 
                for(Vertex currentNode : vertex){ 
                        if(currentNode.name == taxonName){ 
                                taxon = currentNode; 
                        } 
                } 
                ArrayList<Vertex> pathToRoot = new ArrayList<Vertex>(); 
                //Store left or right value from above node... 
                ArrayList<Boolean> isLeft = new ArrayList<Boolean>(); 
                Vertex nodeToAdd = taxon; 
                // create path from taxon to root - first in path is always taxon...last will always be root... 
                pathToRoot.add(taxon); 
                do{ 
                        isLeft.add((nodeToAdd.parent.left==nodeToAdd)?true:false); 
                        nodeToAdd = nodeToAdd.parent; 
                        pathToRoot.add(nodeToAdd); 
                }while(nodeToAdd != root); 
                // Now we have a list of vertices from the taxon to the root and an associate list of the side on which the taxon is of its parent....except for root of course...! 
                 
                // make the non root side have parent of last node in list... 
                // set root as above second vertex 
                //Delete the root 
        } 
         
    } 

}