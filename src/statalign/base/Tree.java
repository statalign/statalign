package statalign.base;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import statalign.base.hmm.Hmm2;
import statalign.base.hmm.HmmNonParam;
import statalign.base.hmm.HmmSilent;
import statalign.base.hmm.HmmTkf92;
import statalign.base.thread.Stoppable;
import statalign.io.DataType;
import statalign.model.score.SubstitutionScore;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Dayhoff;
import statalign.postprocess.plugins.contree.CNetwork;

/**
 * This is the current tree in the MCMC run.
 * It extends Stoppable, so the MCMC run can be terminated in graphical interface mode
 * while it calls a function in the Tree class.
 * @author miklos, novak
 */
public class Tree extends Stoppable implements DataType {

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

    /** The heat parameter for this MCMC chain. */
	public double heat = 1.0d;

    public CNetwork network;

    //boolean changingTree = false;
    boolean changingTree = true;
    
    public String[] previousFullAlign;

    public Tree() {
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
                column.upp = new double[substitutionModel.e.length];
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
            root.calcFelsenRecursively();
            root.calcIndelLogLikeRecursively();
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

    void markAllVerticesUnselected() {
    	for (Vertex v : vertex) {
    		v.selected = false;
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
    
    public void checkPointers() {
    	for (int i = 0; i < vertex.length; i++) {
			if (vertex[i].left != null && vertex[i].right != null) {				
				//System.out.println("Checking vertex "+i);
				vertex[i].checkPointers();				
				AlignColumn p;
				// checking pointer integrity
				for (AlignColumn c = vertex[i].left.first; c != null; c = c.next) {
					p = vertex[i].first;
					while (c.parent != p && p != null)
						p = p.next;
					if (p == null)
						throw new Error(
								"child does not have a parent!!!"
										+ vertex[i] + " "
										+ vertex[i].print(4));
				}
				for (AlignColumn c = vertex[i].right.first; c != null; c = c.next) {
					p = vertex[i].first;
					while (c.parent != p && p != null)
						p = p.next;
					if (p == null)
						throw new Error(
								"child does not have a parent!!!"
										+ vertex[i] + " "
										+ vertex[i].print(4));
				}

			}
		}
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

    public double getLogLike() {
        return root.indelLogLike + root.orphanLogLike;
    }
    /**
     * Recomputes the likelihood from scratch and checks whether everything was correct before.
     * Throws error if it finds inconsistency.
     */
    public void recomputeCheckLogLike() {
    	root.recomputeCheckLogLike();
     }
    
    public int countLeaves() {
        return root.countLeaves();
    }
    
    void checkUppFelsProducts() {
    	double[] result = new double[vertex.length];
    	System.out.println("upp · fels = ");
VERTEX:	for (Vertex v : vertex) {
    		boolean isRoot = (v == root);
    		result[v.index] = 0.0;
    		AlignColumn p = null, pp = null; 
    		if (!isRoot) p = v.parent.first;
    		AlignColumn c = v.first;
    		double tmp = 0;
COLUMN:    	while (c != v.last) {
				
				tmp = Math.log(Utils.calcEmProb(c.seq,c.upp));		
				if (v.left!=null && v.right!=null) {
					//System.out.println(c.seq[0]+" "+c.left.seq[0]+" "+v.left.charTransMatrix[0][0]);
					boolean match = Utils.calcFelsenWithCheck(c.seq,
							c.left!=null?c.left.seq:null,v.left.charTransMatrix,
							c.right!=null?c.right.seq:null,v.right.charTransMatrix);
					if (!match) {
						throw new RuntimeException("Didn't match.");
					}
				}
				result[v.index] += tmp;
				char par = '-', car = ' ';
				if (c != v.last) car = c.mostLikely();
				if (p != null && p != v.parent.last) par = p.mostLikely(); else par = '-';
				System.out.print("\t"+tmp+" ("+par+" -> "+car+") ");
				for (int i=0; i<c.seq.length; i++) {
					System.out.print(c.seq[i]+"/"+c.upp[i]+" ");
				}
				System.out.println();
				
				if (isRoot) { c = c.next; continue COLUMN; }
				
				p = c.parent;
				if (p == null) throw new RuntimeException ("Start p is null.");
				c = c.next; 	
    			while (c.orphan) { 
    				if (c != v.last) car = c.mostLikely(); else car = ' ';
    				tmp = Math.log(Utils.calcEmProb(c.seq,c.upp));
    				result[v.index] += tmp; 
    				System.out.print("\t"+tmp+" ("+"- -> "+car+") ");
	    			for (int i=0; i<c.seq.length; i++) {
	    				System.out.print(c.seq[i]+"/"+c.upp[i]+" ");
	    			}
	    			System.out.println();
    				c = c.next;
    			}
    			pp = c.parent;
    			p = p.next;
    			//if (p != null && p != v.parent.last) System.out.print(" ("+p.mostLikely()+","+pp.mostLikely()+") ");
	    		// Include contributions from all parent columns that 
	    		// experience a deletion on this branch.		    	
	    		while (p != pp) {		 
    				if (p != null && p != v.parent.last) par = p.mostLikely();
    				
    				tmp = Math.log(Utils.calcEmProb(p.seq, p.upp));
	    			result[v.index] += tmp;
	    			System.out.print("\t"+tmp+" ("+par+"-> -) ");		
	    			for (int i=0; i<p.seq.length; i++) {
	    				System.out.print(p.seq[i]+"/"+p.upp[i]+" ");
	    			}
	    			System.out.println();
	    			p = p.next;
	    			if (p == v.parent.last) break COLUMN;
	    		}
    			    		
	    	}
	    	System.out.println("\t"+result[v.index]+" \n");
    	}
    	System.out.println();
    }


    public int getTopVertexId(int ind) {
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
    //public State getState() {
    //	return getState(root);
    //}
	//public State getState(Vertex subtreeRoot) {
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
    	StringBuffer b = new StringBuffer();
        root.print(b);
        return b.toString();
    }
    
	public double getLogPrior() {		
		return 0;
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
    
    @Override
    public boolean perSequenceData() {
    	return false;
    }
    
    @Override
    public String getSummaryAssociatedWith(String sequenceName) {
    	return null;
    }
    
    @Override
    public void removeDataAssociatedWith(String sequenceName) {
    }

//    /** This function is only for testing purposes. */
//    public static void main(String[] args) {
//        /*
//      Tree tree = new Tree(new String[] {"kkkkkk", "kkkkkkkk", "kkkkkkkk", "kkkkkkkk"},
//                   new String[] {"A","B","C", "D"},
//                   new NonParametric(),
//                   new Blosum62());
//      System.out.println(tree.printedTree());
//      tree.root.calcFelsRecursively();
//      System.out.println(tree.root.orphanLogLike);
//           */
//        try {
//            Tree tree = new Tree(new File("newick"), new Blosum62());
//            System.out.println(tree.printedTree());
//        } catch (FileNotFoundException e) {
//        } catch (IOException e) {
//        }
//    }
    
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