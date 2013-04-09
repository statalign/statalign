package statalign.base;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.FieldPosition;
import java.util.ArrayList;

import statalign.ui.ErrorMessage;

/**
 * This is a vertex of the tree.
 * The 'hardcore' functions are implemented in this class, developers are suggested
 * not change functions in it. The implemented functions are quite unreadable, since we
 * opted for efficiency and not readability. You should be able to develop novel
 * functionality of the software package (postprocessing, substitution models, etc.)
 * without touching this class.
 * @author miklos, novak
 */
public class Vertex {

    static final int ROUNDING = 100; // tells the number of digits in rounding when printing the tree
    static final double SELECTING = 0.5; /*this is the probability for selecting a 
    										homologous column not to be changed during topology changing
	 */
    static final double PROPOSAL_HEAT = 1.0;
    static final double EMPTY_WINDOW = 0.01; /* this is the probability that an empty window will be realigned*/
    
    Tree owner;

    Vertex old;

    /**
     * The name of the sequence associated to the vertex.
     * Used only if the vertex is a leaf of the tree.
     */
    public String name;
    
    /** Index of this vertex in the Tree.vertex array */
    public int index;

    /** This reference points to the parent of the vertex */
    public Vertex parent;
    /**
     * This reference points to the left child of the vertex. If the
     * vertex is a leaf, it is set to null.
     */
    public Vertex left;
    /**
     * This reference points to the right child of the vertex. If the
     * vertex is a leaf, it is set to null.
     */
    public Vertex right;

    int length;					// sequence length
    public AlignColumn first;			// first alignment column of Vertex
    AlignColumn last;			// last, virtual alignment column of Vertex (always present)
    String seq;					// original sequence of this Vertex (given for leaves only)

    int winLength;                    // length of window
    AlignColumn winFirst;        // first alignment column of window
    AlignColumn winLast;        // first alignment column past window end
    public boolean selected;                // shows if vertex is part of the selected subtree

    /** The length of the edge that connects this vertex with its parent. */
    public double edgeLength;                            // length of edge to parent vertex
    double[][] charTransMatrix;            // precalculated character transition likelihoods (subst. model)
    double[][] charPropTransMatrix;        // precalculated character transition likelihoods for proposals (subst. model)
    double[][] hmm2TransMatrix;            // precalculated state transition likelihoods for 2-seq HMM (indel model)
    double[][] hmm2PropTransMatrix;        // precalculated state transition likelihoods for 2-seq HMM used for proposals (indel model)
    double[][] hmm3TransMatrix;            // precalculated state transition likelihoods for 3-seq HMM (indel model)
    double[][] hmm3RedTransMatrix;    // precalculated st. trans. likelihoods for 3-seq HMM, silent st. removed (indel model)

    /**
     * The log-sum of the Felsenstein's likelihoods of characters that are inserted into the
     * sequence of this vertex.
     */
    public double orphanLogLike;        // log-sum of the likelihood of each orphan column in subtree (incl. this vertex)
    /**
     * The log-sum of the cumulative insertion-deletion loglikelihoods up to this vertex (ie. summed over the
     * subtree below this vertex.).
     */
    public double indelLogLike;

    public int leafCount;

    Vertex() {
    }
    
	public Vertex(Tree tree) {
		owner = tree;
		old = new Vertex();
	}

	public Vertex(Vertex parent) {
		this.parent = parent;
		owner = parent.owner;
		old = new Vertex();
	}

    Vertex(Tree owner, double edgeLength) {
        this.owner = owner;
        this.edgeLength = edgeLength;
        indelLogLike = 0.0;
        old = new Vertex();
    }

    Vertex(Tree owner, int index, double edgeLength, String name) {
        this.owner = owner;
        this.index = index;
        this.edgeLength = edgeLength;
        this.name = name;
        old = new Vertex();
    }
    
    void addSeq(int[][] seq, String origSeq) {
    	this.seq = origSeq;
        int size = owner.substitutionModel.e.length;
        first = new AlignColumn(this);
        first.seq = new double[size];
        for (int i = 0; i < first.seq.length; i++) {
            first.seq[i] = seq[0][i];
        }
        first.prev = null;
        AlignColumn prev = first;
        for (int i = 1; i < seq.length; i++) {
            AlignColumn actual = new AlignColumn(this);
            prev.next = actual;
            actual.prev = prev;
            actual.seq = new double[size];
            for (int j = 0; j < first.seq.length; j++) {
                actual.seq[j] = seq[i][j];
            }
            prev = actual;
        }
        last = new AlignColumn(this);
        last.prev = prev;
        prev.next = last;
        length = seq.length;
        edgeChangeUpdate();
        indelLogLike = 0.0;
    }
    
    /**
     * Adds empty sequence to an ancestor to allow 3-way alignment.
     * Left and right vertices must exist.
     */
    void addEmptySeq() {
        AlignColumn fake = new AlignColumn(this);
        first = fake;
        last = fake;
        length = 0;
        fake.left = left.last;
        fake.right = right.last;
        left.last.parent = fake;
        left.last.orphan = false;
        right.last.parent = fake;
        right.last.orphan = false;
    }


    Vertex(String descriptor, Vertex parent) {

        /* descriptor is a string describing the nodes below this node */

        //System.out.println(descriptor);

        this.parent = parent;
        owner = parent.owner;
        old = new Vertex();

        /*
	boolean containsParenthesis = false;
	for(int i = 0; i < descriptor.length(); i++){
	    if(descriptor.charAt(i) == '(' || descriptor.charAt(i) == ')'){
		containsParenthesis = true;
		break;
	    }
	}
		 */

        if (descriptor.charAt(0) == '(') {
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
            Vertex leftChild = new Vertex(leftDescriptor, this);
            left = leftChild;

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
            Vertex rightChild = new Vertex(rightDescriptor, this);
            right = rightChild;
            String tempStringlength = "";
            for (j += 1; j < descriptor.length(); j++) {
                tempStringlength += descriptor.charAt(j);
            }
            //System.out.println(tempStringlength);
            if (tempStringlength.length() > 0) {
                edgeLength = Math.max(Double.parseDouble(tempStringlength), 0.01);
            } else {
                edgeLength = 0.01;
            }
            //  edgeLength = 0.01;
            //else we have the root node, which does not need an edgeLength
        } else {
            left = null;
            right = null;
            name = "";
            int i = 0;
            while (descriptor.charAt(i) != ':') {
                name += descriptor.charAt(i);
                i++;
            }
            String tempString = "";
            for (i++; i < descriptor.length(); i++) {
                tempString += descriptor.charAt(i);
            }
            //System.out.println(tempString);
            edgeLength = Math.max(Double.parseDouble(tempString), 0.01);

        }

        // setting the transition matrices
        //	edgeChangeUpdate();

    }

    /**
     * @return the "brother" of this node, i.e. the other descendant of its parent
     */
    public Vertex brother() {
        return parent.left == this ? parent.right : parent.left;
    }

    public void updateTransitionMatrix() {
        //	System.out.println("owner: "+owner);
        //System.out.println("");
        charTransMatrix = owner.substitutionModel.updateTransitionMatrix(charTransMatrix, edgeLength);
        charPropTransMatrix = new double[charTransMatrix.length][];


//        double heat = owner.heat;
//        for (int i = 0; i < charTransMatrix.length; i++) {
//            charPropTransMatrix[i] = charTransMatrix[i].clone();
//			double tempSum = 0;
//			for (int j = 0; j < charPropTransMatrix[i].length; j++) {
//				charPropTransMatrix[i][j] = Math.pow(charPropTransMatrix[i][j], heat);
//				tempSum += charPropTransMatrix[i][j];
//			}
//			for (int j = 0; j < charPropTransMatrix[i].length; j++) {
//
//				charPropTransMatrix[i][j] /= tempSum;
//			}
//
//        }

    }

    void updateHmm2Matrix() {
        hmm2TransMatrix = owner.hmm2.preCalcTransMatrix(hmm2TransMatrix, edgeLength);

        // This is commented out for now. Uncomment to use the rescaling of the transition matrices.
        // Do the same for {@link #updateHmm3Matrix()}.
        
//        double heat = owner.heat;
        double heat = PROPOSAL_HEAT;
        hmm2PropTransMatrix = new double[hmm2TransMatrix.length][];
        for (int i = 0; i < hmm2TransMatrix.length; i++) {
            hmm2PropTransMatrix[i] = hmm2TransMatrix[i].clone();
            /**/
			double tempSum = Utils.log0;
			for (int j = 0; j < hmm2PropTransMatrix[i].length; j++) {

				hmm2PropTransMatrix[i][j] = hmm2PropTransMatrix[i][j]* heat;

				tempSum = Utils.logAdd(tempSum, hmm2TransMatrix[i][j]);
				
			}
			if(tempSum != Double.NEGATIVE_INFINITY){
				for (int j = 0; j < hmm2PropTransMatrix[i].length; j++) {

					hmm2PropTransMatrix[i][j] = hmm2PropTransMatrix[i][j]-tempSum;
				}
			}
            /**/
        }


    }

    void updateHmm3Matrix() {
        if (left != null && right != null) {
            hmm3TransMatrix = owner.hmm3.preCalcTransMatrix(hmm3TransMatrix, left.edgeLength, right.edgeLength);
            hmm3RedTransMatrix = owner.hmm3.preCalcRedTransMatrix(hmm3RedTransMatrix, hmm3TransMatrix);

            
			double heat = PROPOSAL_HEAT;
			
//			double heat = owner.heat;
//			
			for (int i = 0; i < hmm3RedTransMatrix.length; i++) {
				double tempSum = Utils.log0;
				for (int j = 0; j < hmm3RedTransMatrix[i].length; j++) {
					hmm3RedTransMatrix[i][j] = hmm3RedTransMatrix[i][j] * heat;
					tempSum = Utils.logAdd(tempSum, hmm3RedTransMatrix[i][j]);

				}
				if(tempSum != Double.NEGATIVE_INFINITY){
					for (int j = 0; j < hmm3RedTransMatrix[i].length; j++) {
						hmm3RedTransMatrix[i][j] = hmm3RedTransMatrix[i][j] - tempSum;
					}
				}
			}
            
        }
    }

    public void updateHmmMatrices() {
        if (parent != null)
            updateHmm2Matrix();
        updateHmm3Matrix();
    }

    public void setEdgeLength(double x) {
		edgeLength = x;
	}
    
	public static void printChildren(Vertex v){
		if(v.left != null){
			System.out.println(v.index + " " + v.left.index);
			System.out.println(v.index + " " + v.right.index);
			printChildren(v.left);
			printChildren(v.right);
		}
	}

	public static void printEdges(Vertex v){
		System.out.println(v.index + " " + v.edgeLength);
		if(v.left != null){
			printEdges(v.left);
			printEdges(v.right);
		}
	}
    
    public void edgeChangeUpdate() {
        if (parent != null) {
            updateTransitionMatrix();
            updateHmm2Matrix();
            parent.updateHmm3Matrix();
        }
    }

    void fullWin() {
        winFirst = first;
        winLast = last;
        winLength = length;
    }

    static DecimalFormat printDf;
    static {
    	DecimalFormatSymbols dfs = new DecimalFormatSymbols(); 
    	dfs.setDecimalSeparator('.');
    	printDf = new DecimalFormat("0.#####", dfs);
    }

    public void print(StringBuffer b) {

        // print this node and nodes below in bracket notation

        if (left == null && right == null) {
        	b.append(name.replaceAll(" ", ""));
        	b.append(':');
        	printDf.format(edgeLength, b, new FieldPosition(1));
        } else {
        	b.append('(');
        	left.print(b);
        	b.append(',');
        	right.print(b);
        	b.append(')');
            if (parent != null) {
            	b.append(':');
            	printDf.format(edgeLength, b, new FieldPosition(1));
            } else {
                b.append(';');
            }
        }
    }

    public String print(int digits) {

        String pattern = "0.";
        for (int i = 0; i < digits; i++) {
            pattern += "#";
        }
        // print this node and nodes below in bracket notation
        DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator('.');
        DecimalFormat df = new DecimalFormat(pattern, dfs);

        if (left == null && right == null) {
            String x = df.format(edgeLength, new StringBuffer(), new FieldPosition(1)).toString();

            return name.replaceAll(" ", "") + (digits == 0 ? "" : ":" + x);
        } else {
            if (parent != null) {
                //	String x = Double.toString(edgeLength);
                String x = df.format(edgeLength, new StringBuffer(), new FieldPosition(1)).toString();
                return "(" + left.print(digits) + "," + right.print(digits) + ")" + (digits == 0 ? "" : ":" + x);
            } else {
                return "(" + left.print(digits) + "," + right.print(digits) + ");";
            }
        }
    }

    /**
     * This function calculates the Felsenstein likelihood. The result is stored
     * in the orphanLogLike at the root.
     */
    public void calcFelsRecursively() {
        if (left != null && right != null) {
            // System.out.println("calling the left child");
            left.calcFelsRecursively();
            //System.out.println("calling the right child");
            right.calcFelsRecursively();
        }
        calcFelsen();
        calcOrphan();
    }
    /**
     * This function calculates the upper probability vectors, which contain
     * the partial likelihoods for everything except this subtree.
     */
    public void calcUpperRecursively() {
    	calcUpper();
    	if (left != null && right != null) {
            // System.out.println("calling the left child");
            left.calcUpperRecursively();
            //System.out.println("calling the right child");
            right.calcUpperRecursively();
        }
    }
    void calcUpper() {
    	AlignColumn v; 
    	Vertex brother = null;
    	double[][] trans = null;
    	if (this != owner.root) {
    		brother = brother();
    		trans = brother.charTransMatrix;
    	}
    	
        for (v = first; v != last; v = v.next) {
        	if (this == owner.root || v.orphan) {
        		v.upp = owner.substitutionModel.e.clone();
        	}
        	else {
        		AlignColumn p = v.parent;
        		AlignColumn b = (this == parent.left) ? p.right : p.left;
        		v.upp = p.upp.clone();
//        		System.out.print("v ");
//        		for (int i=0; i<v.upp.length; i++) {
//    	          	System.out.print(v.upp[i]+" ");
//    	        }
//    	        System.out.println();
        		if (b != null) {
	        		for (int i=0; i<v.upp.length; i++) {
	        			for (int j=0; j<v.upp.length; j++) {
	        				v.upp[i] += b.seq[j] * trans[i][j];
	        			}
	        		}
//	        		System.out.print("b ");
//	        		for (int i=0; i<v.upp.length; i++) {
//	    	          	System.out.print(b.seq[i]+" ");
//	    	        }
//	    	        System.out.println();
//	        		System.out.print("bv ");
//	        		for (int i=0; i<v.upp.length; i++) {
//	    	          	System.out.print(v.upp[i]+" ");
//	    	        }
//	    	        System.out.println();
        		}
        	}
//    	   for (int i=0; i<v.upp.length; i++) {
//           	System.out.print(v.upp[i]+" ");
//           }
//           System.out.println();
        }
     
    }


    void recomputeCheckLogLike() {
        if (left != null && right != null) {
            // System.out.println("calling the left child");
            left.recomputeCheckLogLike();
            //System.out.println("calling the right child");
            right.recomputeCheckLogLike();
        }
        calcFelsenWithCheck();
        calcOrphanWithCheck();
        calcIndelLogLike(true);
    }


    /** Calculates Felsenstein likelihoods of `this' */
    void calcFelsen() {
        if (left != null && right != null) {
            AlignColumn p;
            AlignColumn l = left.first;
            AlignColumn r = right.first;
            double[] fel1, fel2;
            for (p = first; p != last; p = p.next) {
                fel1 = null;
                fel2 = null;
                while (l != left.last && l.orphan)
                    l = l.next;
                while (r != right.last && r.orphan)
                    r = r.next;
                if (l.parent == p) {
                    fel1 = l.seq;
                    l = l.next;
                }
                if (r.parent == p) {
                    fel2 = r.seq;
                    r = r.next;
                }
                Utils.calcFelsen(p.seq, fel1, left.charTransMatrix, fel2, right.charTransMatrix);
            }
        }
    }

    /** Calculates Felsenstein likelihoods of `this' */
    void calcFelsenWithCheck() {
        if (left != null && right != null) {
            AlignColumn p;
            AlignColumn l = left.first;
            AlignColumn r = right.first;
            double[] fel1, fel2;
            for (p = first; p != last; p = p.next) {
                fel1 = null;
                fel2 = null;
                while (l != left.last && l.orphan)
                    l = l.next;
                while (r != right.last && r.orphan)
                    r = r.next;
                if (l.parent == p) {
                    fel1 = l.seq;
                    l = l.next;
                }
                if (r.parent == p) {
                    fel2 = r.seq;
                    r = r.next;
                }
                boolean match = Utils.calcFelsenWithCheck(p.seq, fel1, left.charTransMatrix, fel2, right.charTransMatrix);
                if (!match) {
                    new ErrorMessage(null, "Felsenstein does not match! ", true);
                }
            }
        }
    }

    /**
     * Calculates the sum of orphan likelihoods in the (inclusive) subtree of `this'
     * which is then stored in `this.orphanLogLikge' (this implies that alignment to
     * parent must exists at the time of calling)
     * Saves previous orphan likelihood into `old' so it shouldn't be called twice in a row.
     * (But does not save previous Felsensteins of `this' in AlignColumn's `seq'.)
     */
    void calcOrphan() {
        old.orphanLogLike = orphanLogLike;

        //orphan likelihoods
        orphanLogLike = 0.0;

        for (AlignColumn actual = first; actual != last; actual = actual.next) {
            if (actual.parent == null || actual.orphan) {
                orphanLogLike += Math.log(Utils.calcEmProb(actual.seq, owner.substitutionModel.e));
            }
        }

        if (left != null && right != null)
            orphanLogLike += left.orphanLogLike + right.orphanLogLike;
    }

    /**
     * Calculates the sum of orphan likelihoods in the (inclusive) subtree of `this'
     * which is then stored in `this.orphanLogLike' (this implies that alignment to
     * parent must exists at the time of calling)
     * Saves previous orphan likelihood into `old' so it shouldn't be called twice in a row.
     * (But does not save previous Felsensteins of `this' in AlignColumn's `seq'.)
     */
    void calcOrphanWithCheck() {
        double oldorphanLogLike = orphanLogLike;

        //orphan likelihoods
        orphanLogLike = 0.0;

        for (AlignColumn actual = first; actual != last; actual = actual.next) {
            if (actual.parent == null || actual.orphan) {
                orphanLogLike += Math.log(Utils.calcEmProb(actual.seq, owner.substitutionModel.e));
            }
        }

        if (left != null && right != null)
            orphanLogLike += left.orphanLogLike + right.orphanLogLike;

        if (oldorphanLogLike - orphanLogLike > 1e-6) {
            new ErrorMessage(null, "Problem with orphan loglike: fast: " + oldorphanLogLike + " slow: " + orphanLogLike, true);
            int index = 0;
            while (owner.vertex[index] != this) {
                index++;
            }
            //System.out.println("We are in vertex "+index+(name == null ? print(0) : "("+name+")"));
            //System.out.print("Selected vertices: ");
            //for(int i = 0; i < owner.vertex.length; i++){
            //	if(owner.vertex[i].selected){
            //		System.out.print(i+(owner.vertex[i].name == null ? owner.vertex[i].print(0) : "("+owner.vertex[i].name+")")+" ");
            //	}
            //}
            //System.out.println();
        }
    }

    void calcRecursivlyOrphanWithCheck() {
        if (left != null && right != null) {
            left.calcRecursivlyOrphanWithCheck();
            right.calcRecursivlyOrphanWithCheck();
        }
        calcOrphanWithCheck();
    }


    public void calcIndelLikeRecursively() {
        if (left != null && right != null) {
            left.calcIndelLikeRecursively();
            right.calcIndelLikeRecursively();
        }
        calcIndelLogLike(false);
    }

    /**
     * Function to calculate the indel logLikelihood of the subtree of `this' including
     * the indel logLikelihood of the alignment between `this.left'&`this' and `this.right'&`this'.
     * Saves previous logLikelihood into `old' Vertex, so it shouldn't be called twice in a row.
     * Assumes it has been called previously on both `left' and `right'.
     * Result is stored in `this'.
     * @param withCheck if true likelihood is re-calculated and checked against the already stored value
     */
    void calcIndelLogLike(boolean withCheck) {
    	if(!withCheck)
    		old.indelLogLike = indelLogLike;
        double newIndelLogLike = 0.0;
        if (left != null && right != null) {
            newIndelLogLike += left.calcIndelLogLikeUp();
            newIndelLogLike += right.calcIndelLogLikeUp();
        }
        if(!withCheck)
        	indelLogLike = newIndelLogLike;
        else if(Math.abs(indelLogLike - newIndelLogLike) > 1e-5)
        	throw new Error("likelihood inconsistency in calcIndelLogLike - stored likelihood: "+indelLogLike+" recomputed: "+newIndelLogLike);
    }

    /**
     * Function to calculate the indel log-likelihood of the subtree of `this' plus
     * the indel log-likelihood of the alignment between `this' & `this.parent'
     * Assumes log-likelihood of subtree is already pre-calculated and stored in `this'
     * @return counted log-likelihood
     */
    double calcIndelLogLikeUp() {
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        final int emitPatt2State[] = owner.hmm2.getEmitPatt2State();

        double indelLogLikeUp = indelLogLike;

        //System.out.println("--------------------------------------------------");
        //printPointers();

        if (parent != null) {
            AlignColumn c = first, p = parent.first;
            int prevk = START, k;

            while (c != last || p != parent.last) {
                if (c.parent != p) {                            // deletion (* -), pattern code 2
                    k = emitPatt2State[2];
                    p = p.next;
                } else if (c.orphan) {                        // insertion (- *), pattern code 1
                    k = emitPatt2State[1];
                    c = c.next;
                } else {                                                // substitution (* *), pattern code 3
                    k = emitPatt2State[3];
                    p = p.next;
                    c = c.next;
                }
                indelLogLikeUp += hmm2TransMatrix[prevk][k];
                prevk = k;
            }

            indelLogLikeUp += hmm2TransMatrix[prevk][END];
        }

        return indelLogLikeUp;
    }

    private double[][][] hmm3ProbMatrix() {
        double[] equDist = owner.substitutionModel.e;
        int hmm3Parent[] = owner.hmm3.getStateEmit()[0];
        int hmm3Left[] = owner.hmm3.getStateEmit()[1];
        int hmm3Right[] = owner.hmm3.getStateEmit()[2];
        int leftLen = left.winLength, rightLen = right.winLength;
        final int START = owner.hmm3.getStart();
        final int END = owner.hmm3.getEnd();

        double probMatrix[][][];                                                                // DP matrix used for 3-seq HMM alignment
        probMatrix = new double[leftLen + 1][rightLen + 1][END];        // don't reserve space for end state

        double emissionProb, felsen[] = new double[equDist.length], tr;
        AlignColumn l = null, r = null;
        int i, j, k, previ, prevj, prevk;

        /* i: left prefix length, j: right prefix length, k: state 
	   l: left child's actual alignment column, r: that of right child */
        for (i = 0; i <= leftLen; i++) {
            for (j = 0; j <= rightLen; j++) {
                probMatrix[i][j][START] = (i == 0 && j == 0) ? 0.0 : Utils.log0;        // update matrix for start state
                for (k = START + 1; k < END; k++) {
                    previ = i - hmm3Left[k];
                    prevj = j - hmm3Right[k];
                    if (previ >= 0 && prevj >= 0) {
                        if (hmm3Parent[k] != 0) {
                            // there's a parent character in alignment column, but still can be an (even double-) deletion
                            Utils.calcFelsen(felsen, hmm3Left[k] != 0 ? l.seq : null, left.charTransMatrix, hmm3Right[k] != 0 ? r.seq : null, right.charTransMatrix);
                            emissionProb = Utils.calcEmProb(felsen, equDist);
                        } else {
                            // no parent, it's an insertion on either of (but never both) branches
                            emissionProb = Utils.calcEmProb(hmm3Left[k] != 0 ? l.seq : r.seq, equDist);
                        }
                        if (previ == 0 && prevj == 0)
                            tr = hmm3RedTransMatrix[START][k];
                        else {
                            for (tr = Utils.log0, prevk = START + 1; prevk < END; prevk++)
                                tr = Utils.logAdd(tr, probMatrix[previ][prevj][prevk] + hmm3RedTransMatrix[prevk][k]);
                        }
                        probMatrix[i][j][k] = Math.log(emissionProb) + tr;
                    } else
                        probMatrix[i][j][k] = Utils.log0;
                }
                r = (j > 0) ? r.next : right.winFirst;
            }
            l = (i > 0) ? l.next : left.winFirst;
        }

        return probMatrix;
    }

    /**
     * Computes backproposal of alignment between this & this.left & this.right, window sizes matter.
     * @return log of backproposal probability
     */
    double hmm3BackProp() {
        int emitPatt2State[] = owner.hmm3.getEmitPatt2State();
        final int START = owner.hmm3.getStart();
        final int END = owner.hmm3.getEnd();
        final int SILENT = owner.hmm3.getSilent();

        AlignColumn l, r, p;
        double retVal = 0.0;
        int k, previ, prevj, prevk;
        int pattCode;                                    // alignment column pattern's binary code

        double probMatrix[][][] = hmm3ProbMatrix();

        /* backproposal calculation */

        previ = prevj = 0;                        // previous prefix lengths (0 at the beginning)
        prevk = START;                                // previous state (start state at the beginning)
        l = left.winFirst;                        // starting alignment columns (window starts)
        r = right.winFirst;
        p = winFirst;
        int silentNum = 0;                        // number of transitions through silent state since last non-silent state

        while (prevk != END) {
            if (l == left.winLast && r == right.winLast && p == winLast)        // go to end state (virtual pattern code: 8)
                pattCode = 8;
            else if (l.parent == p && l.orphan)                    // left insertion (- * -), emission pattern binary code: 2
                pattCode = 2;
            else if (r.parent == p && r.orphan)                    // right insertion (- - *), emission pattern binary code: 1
                pattCode = 1;
            else if (l.parent != p && r.parent != p) {        // double deletion (* - -), silent state, emission pattern binary code: 4
                silentNum++;
                p = p.next;
                continue;
            } else {
                pattCode = 4;                            // there is a parent in all the remaining cases (* ? ?)
                if (l.parent == p)
                    pattCode += 2;                    // could be left/right deletion (* - *)/(* * -) with codes 5/6 or
                if (r.parent == p)
                    pattCode += 1;                    // "full" alignment column (* * *), binary code 7 
            }

            k = emitPatt2State[pattCode];            // real state number for column

            double bp = Utils.log0;
            for (int bpPrevk = START; bpPrevk < END; bpPrevk++)
                bp = Utils.logAdd(bp, probMatrix[previ][prevj][bpPrevk] + hmm3RedTransMatrix[bpPrevk][k]);
            bp = probMatrix[previ][prevj][prevk] - bp;
            if (silentNum == 0)            // non-reducated transition, skipping silent state
                bp += hmm3TransMatrix[prevk][k];
            else                                        // non-reducated transitions, passing through silent state at least once
                bp += hmm3TransMatrix[prevk][SILENT] + hmm3TransMatrix[SILENT][k] + (silentNum - 1) * hmm3TransMatrix[SILENT][SILENT];
            if (bp > 1e-5) {
                //System.out.println("Pozitiv!");
            }

            retVal += bp;

            silentNum = 0;
            prevk = k;
            if ((pattCode & 4) != 0)                // make a parent's step
                p = p.next;
            if ((pattCode & 2) != 0) {            // make a left child's step
                l = l.next;
                previ++;
            }
            if ((pattCode & 1) != 0) {            // make a right child's step
                r = r.next;
                prevj++;
            }
        }

        return retVal;
    }

    private void saveWin() {
        old.winFirst = winFirst;
        old.winLast = winLast.prev;                // temporary reference to last (real) window element (but could be null)
        old.winLength = winLength;
        old.first = first;
        old.length = length;
        old.orphanLogLike = orphanLogLike;
        old.indelLogLike = indelLogLike;
    }

    /**
     * Toggles parent alignment pointers of `this' between parent's old and new sequence (use `toNew' to select).
     * When called, parent's AlignColumn pointers must always point to the new sequence, regardless of `toNew'.
     */
    void toggleUp(boolean toNew) {
        AlignColumn dp = toNew ? parent.winFirst : parent.old.winFirst;        // "destination" parent
        AlignColumn c = first;
        boolean inWin = false;

        AlignColumn p = parent.first;
        while (p != parent.winLast) {
            if (p == parent.winFirst) {
                if (toNew)
                    p = parent.old.winFirst;
                inWin = true;
            }
            if (c.parent != p) {
                if (inWin)
                    dp = dp.next;
                p = p.next;
            } else {
                if (inWin)
                    c.parent = dp;
                c = c.next;
            }
        }
    }

    /**
     * Samples a new alignment between `this' &amp; `this.left' &amp; `this.right', taking window sizes
     * into account. Creates new ancestor sequence, updates all Vertex-stored suppl. data, <b>does not save</b> old
     * data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    double hmm3Align() {
        int hmm3Parent[] = owner.hmm3.getStateEmit()[0];
        int hmm3Left[] = owner.hmm3.getStateEmit()[1];
        int hmm3Right[] = owner.hmm3.getStateEmit()[2];
        int leftLen = left.winLength, rightLen = right.winLength;
        final int START = owner.hmm3.getStart();
        final int END = owner.hmm3.getEnd();
        final int SILENT = owner.hmm3.getSilent();

        int k, previ, prevj, prevk;
        double probMatrix[][][] = hmm3ProbMatrix();
        MuDouble retVal = new MuDouble(0.0);
        double prJump[] = new double[END];            // no need to have an element for end state

        /* stochastic traceback */

        old.winLength = winLength;
        winLength = 0;
        //	AlignColumn wfp = winFirst.prev;				// must save this in case winFirst == winLast


        previ = leftLen;
        prevj = rightLen;
        AlignColumn l = left.winLast, r = right.winLast;
        AlignColumn p = winLast;

        if (p == null) {
            throw new Error("Starting p is null!!!");
        }

        for (k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = probMatrix[previ][prevj][prevk] + hmm3RedTransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);
            //	System.out.println("prevk: "+prevk);
            //if(prevk == START){
            //  System.out.println("prevk is START("+START+"), left.edgelength: "+left.edgeLength+
            //		       " right.edgeLength "+right.edgeLength);
            //  for(int i = 0; i < prJump.length; i++){
            //	System.out.println(prJump[i]);
            //  }
            //}

            if (Utils.chooseOne(Math.exp(hmm3TransMatrix[prevk][k] - hmm3RedTransMatrix[prevk][k]), retVal) == 0) {
                do {
                    p = new AlignColumn(p, true);
                    winLength++;
                } while (Utils.chooseOne(Math.exp(hmm3TransMatrix[SILENT][SILENT]), retVal) == 1);
            }
            if (hmm3Parent[prevk] != 0) {
                p = new AlignColumn(p, true);
                winLength++;
            }
            if (hmm3Left[prevk] != 0) {
                l = l.prev;
                l.parent = p;
                if (l.parent == null) {
                    throw new Error("l.parent is null!!!");
                }
                l.orphan = hmm3Parent[prevk] == 0;
                if (hmm3Parent[prevk] != 0)
                    p.left = l;
                previ--;
            }
            if (hmm3Right[prevk] != 0) {
                r = r.prev;
                r.parent = p;
                if (r.parent == null) {
                    throw new Error("r.parent is null!!!");
                }
                r.orphan = hmm3Parent[prevk] == 0;
                if (hmm3Parent[prevk] != 0)
                    p.right = r;
                prevj--;
            }
        }
        p.setWinFirst(null);
        length += winLength - old.winLength;

        //		System.out.println("saveLeft: "+saveLeft+" saveRight: "+saveRight+" length: "+length);


        // 	left.calcOrphan();
        // 	right.calcOrphan();
        // 	calcFelsen();
        // 	calcIndelLogLike();
        //	System.out.println("printing the alignments at the root of this subtree: "+print());
        //System.out.println("pointers");
        //String[] s = left.printedAlignment();
        //System.out.println(s[0]+"\n"+s[1]+"\n");
        //s = right.printedAlignment();
        //System.out.println(s[0]+"\n"+s[1]+"\n");

        return retVal.value;
    }


    /**
     * Samples a new alignment between `this' &amp; `this.left' &amp; `this.right', taking window sizes
     * into account. Creates new ancestor sequence, updates all Vertex-stored suppl. data, saving old
     * data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    double hmm3AlignWithSave() {
        int hmm3Parent[] = owner.hmm3.getStateEmit()[0];
        int hmm3Left[] = owner.hmm3.getStateEmit()[1];
        int hmm3Right[] = owner.hmm3.getStateEmit()[2];
        int leftLen = left.winLength, rightLen = right.winLength;
        final int START = owner.hmm3.getStart();
        final int END = owner.hmm3.getEnd();
        final int SILENT = owner.hmm3.getSilent();

        int k, previ, prevj, prevk;
        double probMatrix[][][] = hmm3ProbMatrix();
        MuDouble retVal = new MuDouble(0.0);
        double prJump[] = new double[END];            // no need to have an element for end state

        /* stochastic traceback */

        saveWin();
        winLength = 0;
        AlignColumn wfp = winFirst.prev;                // must save this in case winFirst == winLast

        AlignColumn ol = left.winLast, or = right.winLast;
        boolean saveLeft = false, saveRight = false;
        if (left.left == null || left.right == null || !left.left.selected || !left.right.selected) {
            left.saveWin();
            saveLeft = true;
        }
        if (right.left == null || right.right == null || !right.left.selected || !right.right.selected) {
            right.saveWin();
            saveRight = true;
        }

        previ = leftLen;
        prevj = rightLen;
        AlignColumn l = left.winLast, r = right.winLast;
        AlignColumn p = winLast;

        if (p == null) {
            throw new Error("Starting p is null!!!");
        }

        for (k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = probMatrix[previ][prevj][prevk] + hmm3RedTransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);
            //	System.out.println("prevk: "+prevk);
            //if(prevk == START){
            //  System.out.println("prevk is START("+START+"), left.edgelength: "+left.edgeLength+
            //		       " right.edgeLength "+right.edgeLength);
            //  for(int i = 0; i < prJump.length; i++){
            //	System.out.println(prJump[i]);
            //  }
            //}

            if (Utils.chooseOne(Math.exp(hmm3TransMatrix[prevk][k] - hmm3RedTransMatrix[prevk][k]), retVal) == 0) {
                do {
                    p = new AlignColumn(p, true);
                    winLength++;
                } while (Utils.chooseOne(Math.exp(hmm3TransMatrix[SILENT][SILENT]), retVal) == 1);
            }
            if (hmm3Parent[prevk] != 0) {
                p = new AlignColumn(p, true);
                winLength++;
            }
            if (hmm3Left[prevk] != 0) {
                if (saveLeft) {
                    ol = ol.prev;
                    l = new AlignColumn(l, false);
                    l.saveBoth(ol);
                } else {
                    l = l.prev;
                }
                l.parent = p;
                if (l.parent == null) {
                    throw new Error("l.parent is null!!!");
                }
                l.orphan = hmm3Parent[prevk] == 0;
                if (hmm3Parent[prevk] != 0)
                    p.left = l;
                previ--;
            }
            if (hmm3Right[prevk] != 0) {
                if (saveRight) {
                    or = or.prev;
                    r = new AlignColumn(r, false);
                    r.saveBoth(or);
                } else {
                    r = r.prev;
                }
                r.parent = p;
                if (r.parent == null) {
                    throw new Error("r.parent is null!!!");
                }
                r.orphan = hmm3Parent[prevk] == 0;
                if (hmm3Parent[prevk] != 0)
                    p.right = r;
                prevj--;
            }
        }
        p.setWinFirst(wfp);
        length += winLength - old.winLength;

        //		System.out.println("saveLeft: "+saveLeft+" saveRight: "+saveRight+" length: "+length);

        if (saveLeft) {
            l.setWinFirst(left.winFirst.prev);
            if (left.left != null && left.right != null) {
                left.left.toggleUp(true);
                left.right.toggleUp(true);
            }
        }
        if (saveRight) {
            r.setWinFirst(right.winFirst.prev);
            if (right.left != null && right.right != null) {
                right.left.toggleUp(true);
                right.right.toggleUp(true);
            }
        }

        //	left.printPointers();
        //right.printPointers();

        left.calcOrphan();
        right.calcOrphan();
        calcFelsen();
        calcIndelLogLike(false);
        //		System.out.println("printing the alignments at the root of this subtree: "+print());
        //		System.out.println("pointers");
        //	String[] s = left.printedAlignment();
        //System.out.println(s[0]+"\n"+s[1]+"\n");
        //s = right.printedAlignment();
        //System.out.println(s[0]+"\n"+s[1]+"\n");

        return retVal.value;
    }

    private double[][][] hmm2ProbMatrix() {
        double[] equDist = owner.substitutionModel.e;
        int hmm2Parent[] = owner.hmm2.getStateEmit()[0];
        int hmm2Child[] = owner.hmm2.getStateEmit()[1];
        int parentLen = parent.winLength, childLen = winLength;
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        Vertex brother = brother();

        double probMatrix[][][];                                                                    // DP matrix used for 2-seq HMM alignment
        probMatrix = new double[parentLen + 1][childLen + 1][END];        // don't reserve space for end state

        double emissionProb, felsen[] = new double[equDist.length], tr;
        AlignColumn c = null;                // child
        AlignColumn p = null;                // parent
        AlignColumn b;                            // brother
        int i, j, k, previ, prevj, prevk;

        /* i: parent prefix length, j: child prefix length, k: state  */
        for (i = 0; i <= parentLen; i++) {
//        	if (!owner.changingTree && p != null) {
//	        	double[] partialLikelihood1 = p.upp;
//	        	for (int iii=0; iii<partialLikelihood1.length; iii++) {
//	        		System.out.print(partialLikelihood1[iii]+ " ");
//	        	}
//	        	System.out.println();
//        	}
            for (j = 0; j <= childLen; j++) {
                probMatrix[i][j][START] = (i == 0 && j == 0) ? 0.0 : Utils.log0;        // update matrix for start state
                for (k = START + 1; k < END; k++) {
                    previ = i - hmm2Parent[k];
                    prevj = j - hmm2Child[k];
                    if (previ >= 0 && prevj >= 0) {
                        if (hmm2Parent[k] != 0) {                // parent present: substitution (* *) or deletion (* -)
                            b = (this == parent.left) ? p.right : p.left;
                            Utils.calcFelsen(felsen, hmm2Child[k] != 0 ? c.seq : null, charTransMatrix, b != null ? b.seq : null, brother.charTransMatrix);
                            if (parent == owner.root || p.orphan || owner.changingTree) {
                                emissionProb = Utils.calcEmProb(felsen, equDist);
                            }
                            else {
//                            	for (int ii=0; ii<p.upp.length; ii++) {
//                            		System.out.print(p.upp[ii]+" ");
//                            	}
//                            	System.out.println();
                            	emissionProb = Utils.calcEmProb(felsen, p.upp);
                            }
//                            	double[] partialLikelihood = p.parent.seq;                             	
//                            	for (int ii=0; ii<partialLikelihood.length; ii++) {
//                            		if (partialLikelihood[ii] == 0 || partialLikelihood[ii] == Double.POSITIVE_INFINITY) {
//                            			if (p == null) {
//                            				System.out.println("p is null, don't you know?");
//                            			}
//                            			double[] partialLikelihood1 = p.parent.seq;
//                        	        	for (int iii=0; iii<partialLikelihood1.length; iii++) {
//                        	        		System.out.print(partialLikelihood1[iii]+ " ");
//                        	        	}
//                        	        	System.out.println();
//                        	        	partialLikelihood1 = p.seq;
//                        	        	for (int iii=0; iii<partialLikelihood1.length; iii++) {
//                        	        		System.out.print(partialLikelihood1[iii]+ " ");
//                        	        	}
//                        	        	System.out.println();
//                            			throw new RuntimeException((b==null)+" "+hmm2Child[k]+" "+i+" "+j+" "+k+" "+index+" "+parent.index+"\n");
//                            		}
//                            		if (p.seq[ii] == 0 || p.seq[ii] == Double.POSITIVE_INFINITY) {
//                            			System.out.println("ii = "+ii+", p.seq[ii] = "+p.seq[ii]);
//                            		}
//                            		partialLikelihood[ii] /= p.seq[ii];
//                            	}
//                                emissionProb = Utils.calcEmProb(felsen, partialLikelihood);
//                            }
                        } else {                    // insertion (- *)
                            emissionProb = Utils.calcEmProb(c.seq, equDist);
                        }
                        if (previ == 0 && prevj == 0)
                            tr = hmm2PropTransMatrix[START][k];
                        else {
                            for (tr = Utils.log0, prevk = START + 1; prevk < END; prevk++)
                                tr = Utils.logAdd(tr, probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k]);
                        }
                        probMatrix[i][j][k] = Math.log(emissionProb) + tr;
                    } else
                        probMatrix[i][j][k] = Utils.log0;
                }
                c = (j > 0) ? c.next : winFirst;
            }
            p = (i > 0) ? p.next : parent.winFirst;
        }
        System.out.println();
        return probMatrix;
    }

    /**
     * Computes backproposal of alignment between `this' & `this.parent', window sizes matter.
     * @return log of backproposal probability
     */
    double hmm2BackProp() {
        int emitPatt2State[] = owner.hmm2.getEmitPatt2State();
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();

        AlignColumn c, p;
        double retVal = 0.0;
        int k, previ, prevj, prevk;
        int pattCode;                                    // alignment column pattern's binary code

        double probMatrix[][][] = hmm2ProbMatrix();

        /* backproposal calculation */

        previ = prevj = 0;                        // previous prefix lengths (0 at the beginning)
        prevk = START;                                // previous state (start state at the beginning)
        c = winFirst;                                    // starting alignment columns (window starts)
        p = parent.winFirst;

        while (prevk != END) {
            if (c == winLast && p == parent.winLast)            // go to end state (virtual pattern code 4)
                pattCode = 4;
            else if (c.parent != p)                                             // deletion (* -), pattern code 2
                pattCode = 2;
            else if (c.orphan)                                                        // insertion (- *), pattern code 1
                pattCode = 1;
            else                                                                                // substitution (* *), pattern code 3
                pattCode = 3;

            k = emitPatt2State[pattCode];                                // real state number for column

            double bp = Utils.log0;
            for (int bpPrevk = START; bpPrevk < END; bpPrevk++)
                bp = Utils.logAdd(bp, probMatrix[previ][prevj][bpPrevk] + hmm2PropTransMatrix[bpPrevk][k]);
            bp = probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k] - bp;

            retVal += bp;

            if ((pattCode & 2) != 0) {            // make a parent's step
                p = p.next;
                previ++;
            }
            if ((pattCode & 1) != 0) {            // make a child's step
                c = c.next;
                prevj++;
            }
            prevk = k;
        }

        return retVal;
    }

    /**
     * Samples a new alignment between `this' &amp; `this.parent', taking window sizes into account.
     * <b>Does not update</b> any Vertex-stored suppl. data, <b>does not save</b> old data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    double hmm2Align() {
        int hmm2Parent[] = owner.hmm2.getStateEmit()[0];
        int hmm2Child[] = owner.hmm2.getStateEmit()[1];
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        int parentLen = parent.winLength, childLen = winLength;
        boolean isLeft = parent.left == this;
        //System.out.println("Aligning "+(isLeft ? "left" : "right")+" in hmm2Align");

        int k, previ, prevj, prevk;
        double probMatrix[][][] = hmm2ProbMatrix();
        MuDouble retVal = new MuDouble(0.0);
        double prJump[] = new double[END];                        // no need to have an element for end state

        /* stochastic traceback */

        //	parent.saveWin();
        parent.winLength = 0;

        //	AlignColumn oc = winLast;				       			// original child
        //	AlignColumn op = parent.winLast;		     				// original parent

        previ = parentLen;
        prevj = childLen;
        AlignColumn c = winLast;                                   // child
        AlignColumn p = parent.winLast;                                  // parent

        for (k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);

            if (hmm2Parent[prevk] != 0) {
                p = p.prev;
                /*
		  if(isLeft){
		  p.left = null;
		  }
		  else{
		  p.right = null;
		  }
				 */
                previ--;
                parent.winLength++;
            }
            if (hmm2Child[prevk] != 0) {
                c = c.prev;
                c.parent = p;
                c.orphan = hmm2Parent[prevk] == 0;
                if (hmm2Parent[prevk] != 0) {
                    if (isLeft)
                        p.left = c;
                    else
                        p.right = c;
                }
                prevj--;
            }
        }

        //	p.setWinFirst(parent.winFirst.prev);
        //	assert (parent.winLength == parent.old.winLength) : "HMM2 alignment error";

        //for(c = winFirst.prev; c != null && c.parent == parent.old.winFirst; c = c.prev)
        //  c.parent = parent.winFirst;

        //calcOrphan();
        //parent.calcFelsen();
        //parent.calcOrphan();
        //parent.calcIndelLogLike();

        return retVal.value;
    }

    /**
     * Samples a new alignment between `this' &amp; `this.parent', taking window sizes into account.
     * Updates all Vertex-stored suppl. data, saving old data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    public double hmm2AlignWithSave() {
        int hmm2Parent[] = owner.hmm2.getStateEmit()[0];
        int hmm2Child[] = owner.hmm2.getStateEmit()[1];
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        int parentLen = parent.winLength, childLen = winLength;
        boolean isLeft = parent.left == this;

        int k, previ, prevj, prevk;
        double probMatrix[][][] = hmm2ProbMatrix();
        MuDouble retVal = new MuDouble(0.0);
        double prJump[] = new double[END];                        // no need to have an element for end state

        /* stochastic traceback */

        parent.saveWin();
        parent.winLength = 0;

        AlignColumn oc = winLast;                                            // original child
        AlignColumn op = parent.winLast;                            // original parent
        boolean saveChild = false;
        if (left == null || right == null || !left.selected || !right.selected) {
            saveWin();
            saveChild = true;
        }

        previ = parentLen;
        prevj = childLen;
        AlignColumn c = winLast;                                            // child
        AlignColumn p = parent.winLast;                                // parent

        for (k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);

            if (hmm2Parent[prevk] != 0) {
                op = op.prev;
                p = new AlignColumn(p, true);
                if (isLeft)
                    p.saveRight(op);
                else
                    p.saveLeft(op);
                p.saveParent(op);
                previ--;
                parent.winLength++;
            }
            if (hmm2Child[prevk] != 0) {
                if (saveChild) {
                    oc = oc.prev;
                    c = new AlignColumn(c, false);
                    c.saveBoth(oc);
                } else {
                    c = c.prev;
                }
                c.parent = p;
                c.orphan = hmm2Parent[prevk] == 0;
                if (hmm2Parent[prevk] != 0) {
                    if (isLeft)
                        p.left = c;
                    else
                        p.right = c;
                }
                prevj--;
            }
        }

        p.setWinFirst(parent.winFirst.prev);
        if (isLeft)
            parent.right.toggleUp(true);
        else
            parent.left.toggleUp(true);
        assert (parent.winLength == parent.old.winLength) : "HMM2 alignment error";

        if (saveChild) {
            c.setWinFirst(winFirst.prev);
            if (left != null && right != null) {
                left.toggleUp(true);
                right.toggleUp(true);
            }
        }

        for (c = winFirst.prev; c != null && c.parent == parent.old.winFirst; c = c.prev)
            c.parent = parent.winFirst;

        calcOrphan();
        parent.calcFelsen();
        parent.calcOrphan();
        parent.calcIndelLogLike(false);

        return retVal.value;
    }

    /**
     * This function returns the number of leaves that are below this vertex.
     * @return the number of leaves that are below this vertex.
     */
    int countLeaves() {
        return leafCount = (left == null && right == null) ? 1 : left.countLeaves() + right.countLeaves();
    }

    /** Selects a subtree */
    public void selectSubtree(double[] weights, int level) {
        selected = true;
        // continue below with prescribed probability
        if (left != null && right != null) {
            if (Utils.generator.nextDouble() < weights[level]) {
                left.selectSubtree(weights, level + 1);
                right.selectSubtree(weights, level + 1);
            } else {
                left.selected = false;
                right.selected = false;
            }
        }
    }

    /** Marks `this' as the last selected Vertex of its subtree. */
    void lastSelected() {
        selected = true;
        if (left != null && right != null) {
            left.selected = false;
            right.selected = false;
        }
    }

    /** This function cuts out a window and realigns in the selected subtree. */
    public double selectAndResampleAlignment() {

    	System.out.println(index+" "+name+" "+edgeLength);
    	System.out.println(owner.hmm2.params[0]+" "+owner.hmm2.params[1]+" "+owner.hmm2.params[2]);
    	System.out.println(owner.printedTree());
    	StringBuffer sb = new StringBuffer();
//    	if (parent != null) {
//    		parent.print(sb);
//    	}
//    	else {
    		print(sb);
//    	}
        System.out.println(sb.toString());
        
        // this code below checks pointer integrity...
        //    for(int i = 0; i < owner.vertex.length - 1; i++){
        //	owner.vertex[i].calcIndelLogLikeUp();//
        //}
        // select the beginning and end of the alignment

        // select the beginning and end of the window
        MuDouble p = new MuDouble(1.0);
        winLength = Utils.linearizerWeight(length, p, Utils.WINDOW_MULTIPLIER*Math.sqrt(length));
        //System.out.print("Length: "+length+"\t");
        //System.out.print("Window length: "+winLength+"\t");
        int b = (length - winLength == 0 ? 0 : Utils.generator.nextInt(length - winLength));
        System.out.println("Initial window length\t "+winLength);
        //System.out.print("b: "+b+"\t");
        AlignColumn actualAC = first;
        for (int i = 0; i < b; i++) {
            actualAC = actualAC.next;
        }
        winFirst = actualAC;
        for (int i = 0; i < winLength; i++) {
            actualAC = actualAC.next;
        }
        winLast = actualAC;
        
       // if (winLength == 33 && index == 3) {
        	owner.root.calcFelsRecursively();
        	owner.root.calcUpperRecursively();
        	owner.changingTree = false;
      //  }
        //double bpp = -Math.log(p.value * (length - winLength == 0 ? 1 : length - winLength));
        double bpp = -Math.log(p.value / (length - winLength + 1));

        System.out.println("bpp after window select\t "+bpp);

        // select window down (recursively)
        if (left != null && right != null) {
        	left.selectWindow();
        	right.selectWindow();
        }
        
        selectWindowUp();
	
        printToScreenAlignment(b,b+winLength);
        // compute alignment backproposal
        //bpp += doRecBackprop();
        double recBackprop = doRecBackprop();
        System.out.println("recBackprop\t "+recBackprop);
        bpp += recBackprop;
        if (parent != null) {
            //bpp += hmm2BackProp();
        	double hmm2BackProp = hmm2BackProp();
        	System.out.println("hmm2BackProp\t "+hmm2BackProp);
        	bpp += hmm2BackProp;
        }

        // align the sequences
        double bppProp = doRecAlign();
        System.out.println("bppProp after doRecAlign()\t "+bppProp);
        if (parent != null) {
        	//bppProp += hmm2AlignWithSave();
        	double hmm2AlignProp = hmm2AlignWithSave();
        	System.out.println("hmm2AlignProp\t "+hmm2AlignProp);
            bppProp += hmm2AlignProp;
            parent.calcAllUp();
        } else {
            calcOrphan();
        }

        bpp += bppProp;
        System.out.println("bpp after proposal\t "+bpp);

        if(Utils.DEBUG) {
        	// check proposal - backproposal consistency
        	double bppBack = doRecBackprop();
        	if(parent != null)
        		bppBack += hmm2BackProp();
        	if(Math.abs(bppProp+bppBack) > 1e-5) {
        		System.out.println("Proposal - backproposal inconsistent! Prop: "+bppProp+" Back: "+bppBack);
        	}
        }

        // backproposal probability for cutting out the window
        //bpp += Math.log(Utils.linearizerWeightProb(length, winLength, Utils.WINDOW_MULTIPLIER*Math.sqrt(length)) 
        //		* (length - winLength));
        System.out.println("Final window length\t "+winLength);
        bpp += Math.log(Utils.linearizerWeightProb(length, winLength, Utils.WINDOW_MULTIPLIER*Math.sqrt(length)) 
               		/ (length - winLength + 1));

    	System.out.println("Total bpp\t "+bpp);
        //System.out.print(" Prop: "+bppProp+" it's doublecheck: "+bppBack+" bpp: "+bpp+" ");
    	printToScreenAlignment(b,b+winLength);
    //	if (winLength == 33 && index == 3) {
        	owner.changingTree = true;
    //    }
        return bpp;
    }
    
    public double realignToParent() {
    	return realignToParent(false);
    }
    
    /**
     * Samples a new alignment between the subtree rooted at this vertex and the rest of the sequences.
     * Only the alignment between this vertex and its parent is changed.
     * @param useCurrentWin if true then keeps the current window selection for this node
     * 		(backproposal and proposal for window selection must be computed by the caller)
     * @return the ratio between backproposal and proposal probabilities
     */
    public double realignToParent(boolean useCurrentWin) {
    	if(parent == null)
    		throw new Error("realignToParent was called on the root vertex");
    	
    	double bpp = 0;
    	
        if(!useCurrentWin) {
        	MuDouble p = new MuDouble(1.0);
	        winLength = Utils.linearizerWeight(length, p, Utils.WINDOW_MULTIPLIER*Math.sqrt(length));
	        //System.out.print("Length: "+length+"\t");
	        //System.out.print("Window length: "+winLength+"\t");
	        int b = (length - winLength == 0 ? 0 : Utils.generator.nextInt(length - winLength));
	        //System.out.print("b: "+b+"\t");
	        AlignColumn actualAC = first;
	        for (int i = 0; i < b; i++) {
	            actualAC = actualAC.next;
	        }
	        winFirst = actualAC;
	        for (int i = 0; i < winLength; i++) {
	            actualAC = actualAC.next;
	        }
	        winLast = actualAC;
	        
	        //double bpp = -Math.log(p.value * (length - winLength == 0 ? 1 : length - winLength));
	        bpp -= Math.log(p.value / (length - winLength + 1));
        }

        // nothing is selected for realignment below this vertex
        lastSelected();

        // window is projected to parent
        selectWindowUp();

        // compute alignment backproposal
        bpp += hmm2BackProp();

        // align the sequences
        double bppProp = hmm2AlignWithSave();
        bpp += bppProp;
        parent.calcAllUp();

        if(Utils.DEBUG) {
        	// check proposal - backproposal consistency
        	double bppBack = hmm2BackProp();
        	if(Math.abs(bppProp+bppBack) > 1e-5) {
        	  System.out.println("Proposal - backproposal inconsistent in realignToParent! Prop: "+bppProp+" Back: "+bppBack);
        	}
        }

        // backproposal probability for cutting out the window
        //bpp += Math.log(Utils.linearizerWeightProb(length, winLength, Utils.WINDOW_MULTIPLIER*Math.sqrt(length)) 
        //		* (length - winLength));
        bpp += Math.log(Utils.linearizerWeightProb(length, winLength, Utils.WINDOW_MULTIPLIER*Math.sqrt(length)) 
               		/ (length - winLength + 1));

        // 	System.out.print(" Prop: "+bppProp+" it's doublecheck: "+bppBack+" bpp: "+bpp+" ");

        return bpp;
    }


    void selectWindowUp() {
        if (parent != null) {
            if (winFirst.prev == null) {
                parent.winFirst = parent.first;
            } else {
                parent.winFirst = (winFirst.prev.orphan ? winFirst.prev.parent : winFirst.prev.parent.next);
            }
            parent.winLast = winLast.parent;

            parent.winLength = 0;
            for (AlignColumn actualAC = parent.winFirst; actualAC != parent.winLast; actualAC = actualAC.next)
                parent.winLength++;
        }
    }
    
    double doRecAlign() {
        if (left != null && right != null && left.selected && right.selected) {
            double ret = left.doRecAlign() + right.doRecAlign();
            ret += hmm3AlignWithSave();
            return ret;
        }
        return 0.0;
    }

    /** Computes recursively the log of the backproposal probability of the selected subtree. */
    double doRecBackprop() {
        if (left != null && right != null && left.selected && right.selected) {
            double ret = left.doRecBackprop() + right.doRecBackprop();
            ret += hmm3BackProp();
            return ret;
        }
        return 0.0;
    }

    /**
     * Restores all the changes an alignment resampling on the currently selected subtree has produced.
     * Must be called on selected subtree root.
     */
    public void alignRestore() {
        doRecRestore();
        if (parent != null) {
            if (parent.left == this)
                parent.right.toggleUp(false);
            else
                parent.left.toggleUp(false);
            for (AlignColumn c = winFirst.prev; c != null && c.parent == parent.winFirst; c = c.prev)
                c.parent = parent.old.winFirst;
            parent.doRestore();
            if (parent.parent != null) {
                boolean isLeft = parent.parent.left == parent;
                for (AlignColumn p = parent.winFirst; p != parent.winLast; p = p.next) {
                    if (!p.orphan) {
                        if (isLeft)
                            p.parent.left = p;
                        else
                            p.parent.right = p;
                    }
                }
            }
            parent.calcAllUp();
        }
    }

    /**
     * Restores old sequences in subtree of `this' recursively and up-alignments in non-selected nodes.
     * Assumes `this' is selected.
     */
    void doRecRestore() {
        if (left != null && right != null) {
            if (left.selected && right.selected) {
                left.doRecRestore();
                right.doRecRestore();
            } else {
                left.toggleUp(false);
                right.toggleUp(false);
            }
        }
        doRestore();
    }

    /**
     * Restores old sequence into `this', pointed to by `old.winFirst' & `old.winLast'.
     * Also restores orphan and indel likelihoods.
     */
    void doRestore() {
        winFirst = old.winFirst;
        winLength = old.winLength;
        length = old.length;
        winLast.prev = old.winLast;
        if (winFirst.prev != null)
            winFirst.prev.next = winFirst;
        else
            first = winFirst;
        orphanLogLike = old.orphanLogLike;
        indelLogLike = old.indelLogLike;
    }

    /** This function recursively selects the boundaries of the sequence that is to be resampled. */
    void selectWindow() {
        if (!selected) {
            return;
        }
        AlignColumn p = parent.first;
        AlignColumn c = first;
        while (p != parent.winFirst) {
            if (c.parent != p) {
                p = p.next;
            } else {
                c = c.next;
            }
        }
        winFirst = c;
        winLength = 0;
        //end of the window
        while (p != parent.winLast || (c.parent == p && c.orphan)) {
            //if(p == parent.winLast){
            //		found = true;
            //	}
            if (c.parent != p) {
                p = p.next;
            } else {
                c = c.next;
                winLength++;
            }
        }
        winLast = c;

        if (left != null && right != null) {
            left.selectWindow();
            right.selectWindow();
        }
    }

    /**
     * Swaps {@code this} with its uncle, {@link #hmm3Align()}'ing parent and grandparent,
     * the latter is {@link #hmm2Align()}'ed as well.
     * Assumes {@code this} has a non-null grandparent.
     * @return log-quotient of backproposal and proposal
     */
    double swapWithUncle1() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;
        double ret = 0.0;

        fullWin();
        brother().fullWin();
        uncle.fullWin();
        parent.fullWin();
        grandpa.fullWin();
        if (grandpa.parent != null)
            grandpa.parent.fullWin();

        lastSelected();
        brother().lastSelected();
        uncle.lastSelected();
        parent.selected = true;
        grandpa.selected = true;

        ret += parent.hmm3BackProp();
        ret += grandpa.hmm3BackProp();
        if (grandpa.parent != null)
            ret += grandpa.hmm2BackProp();

        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;

        ret += uncle.parent.hmm3AlignWithSave();
        ret += grandpa.hmm3AlignWithSave();
        if (grandpa.parent != null) {
            ret += grandpa.hmm2AlignWithSave();
            grandpa.parent.calcAllUp();
        } else {
            grandpa.calcOrphan();
        }

        return ret;
    }

    /**
     * Restores the exact state just before the call of {@link #swapWithUncle1()}.
     * Must be called on ex-uncle node.
     * Assumes {@code this} has a non-null grandparent.
     */
    void swapBackUncle1() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;

        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;

        grandpa.alignRestore();
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Swaps {@code this} with its uncle, but proposes new alignments only between this node and its parent, and
     * the uncle node and its parent, every other alignment is kept fixed. Slow because full sequence
     * alignment is done.
     * Assumes {@code this} has a non-null grandparent.
     * @return log-quotient of backproposal and proposal
     */
    public double swapWithUncleAlignToParent() {
    	owner.changingTree = true;
        Vertex uncle = parent.brother(), grandpa = parent.parent;
        double ret = 0.0;

        // make node and window selection
        lastSelected();
        uncle.lastSelected();
        
        fullWin();
        parent.fullWin();
        uncle.fullWin();
        grandpa.fullWin();
        
        // compute alignment backproposal
        ret += hmm2BackProp();
        ret += uncle.hmm2BackProp();

        System.out.println("log back proposal probability = "+ret);
        // do the swap
        parentNewChild(uncle);			// order is important here
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;

        // align the sequences (uncle first, that's now below!)
        double bppProp = uncle.hmm2AlignWithSave();
        bppProp += hmm2AlignWithSave();
        System.out.println("log forward proposal probability = "+(-bppProp));
        ret += bppProp;
        parent.calcAllUp();

        if(Utils.DEBUG) {
        	// check proposal - backproposal consistency
        	double bppBack = uncle.hmm2BackProp();
        	bppBack += hmm2BackProp();
        	if(Math.abs(bppProp+bppBack) > 1e-5) {
        	  System.out.println("Proposal - backproposal inconsistent in swapWithUncleAlignToParent! Prop: "+bppProp+" Back: "+bppBack);
        	}
        }
        //owner.changingTree = false;
        return ret;
	}

    /**
     * Restores the exact state just before the call to {@link #swapWithUncleAlignToParent()}.
     * Must be called on ex-uncle node.
     * Assumes {@code this} has a non-null grandparent.
     */
   public void swapBackUncleAlignToParent() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;

        // swap back
        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;

        // restore alignment
        alignRestore();
        uncle.alignRestore();
    }

	/**
     * Swaps {@code this} with its uncle, but keeps alignment columns where {@code this}, its parent,
     * the uncle node and its parent are all aligned unchanged with a given probability, and only
     * proposes new alignments between these selected alignment anchors columns. Alignments between
     * other sequences than the above mentioned four are kept unchanged.
     * Assumes {@code this} has a non-null grandparent.
     * @return log-quotient of backproposal and proposal
     */
    public double fastSwapWithUncle() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;
        double ret = 0.0;
        Vertex starter;

        //System.out.println("fast swap here"+grandpa.print());

        lastSelected();
        brother().lastSelected();
        uncle.lastSelected();
        parent.selected = true;
        grandpa.selected = true;

        setAllAlignColumnsUnselected();
        brother().setAllAlignColumnsUnselected();
        uncle.setAllAlignColumnsUnselected();
        parent.setAllAlignColumnsUnselected();
        grandpa.setAllAlignColumnsUnselected();

        if (grandpa.parent != null) {
            grandpa.parent.selected = true;
            grandpa.brother().selected = false;
            grandpa.parent.setAllAlignColumnsUnselected();
            starter = grandpa.parent;
            //  System.out.println("has grandgrandpa");
        } else {
            starter = grandpa;
            //System.out.println("does not have grandgrandpa");
        }

        ret += starter.selectAnchors();
        //System.out.println("RET after selecting the anchors: "+ret);
        //System.out.println("Backproposing this would be:    "+starter.backproposeAnchors());

        //backproposals
        //for parent
        AlignColumn p = parent.last, pf = parent.last;
        AlignColumn t = last, tf = last;
        AlignColumn b = brother().last, bf = brother().last;
        while (p != null) {
            // System.out.println(".");
            pf = pf.findWindowStart();
            tf = tf.findWindowStart();
            bf = bf.findWindowStart();
            if (p.prev != pf || p.emptyWindow || true) {
                parent.winFirst = pf;
                parent.winLast = p;
                winFirst = tf;
                winLast = t;
                brother().winFirst = bf;
                brother().winLast = b;
                ret += parent.hmm3BackProp();
            }
            p = pf.prev;
            pf = pf.prev;
            t = tf.prev;
            tf = tf.prev;
            b = bf.prev;
            bf = bf.prev;
        }
        //for grandpa
        AlignColumn g = grandpa.last, gf = grandpa.last;
        p = parent.last;
        pf = parent.last;
        AlignColumn u = uncle.last, uf = uncle.last;
        while (g != null) {
//			  System.out.println(".");
            gf = gf.findWindowStart();
            pf = pf.findWindowStart();
            uf = uf.findWindowStart();
            if (g.prev != gf || g.emptyWindow || true) {
                grandpa.winFirst = gf;
                grandpa.winLast = g;
                parent.winFirst = pf;
                parent.winLast = p;
                uncle.winFirst = uf;
                uncle.winLast = u;
                ret += grandpa.hmm3BackProp();
            }
            g = gf.prev;
            gf = gf.prev;
            p = pf.prev;
            pf = pf.prev;
            u = uf.prev;
            uf = uf.prev;

        }
        // if there is a grand-grandpa...
        if (grandpa.parent != null) {
            g = grandpa.last;
            gf = grandpa.last;
            AlignColumn gg = grandpa.parent.last, ggf = grandpa.parent.last;
            while (gg != null) {
                ggf = ggf.findWindowStart();
                gf = gf.findWindowStart();
                if (gg.prev != ggf || gg.emptyWindow || true) {
                    grandpa.parent.winFirst = ggf;
                    grandpa.parent.winLast = gg;
                    grandpa.winFirst = gf;
                    grandpa.winLast = g;
                    ret += grandpa.hmm2BackProp();
                }
                g = gf.prev;
                gf = gf.prev;
                gg = ggf.prev;
                ggf = ggf.prev;
            }
        }

        //System.out.println("RET after alignment backproposal: "+ret);

        // 	/////////////////////////////////////////////////
        // 	System.out.println("Printing the pointers +alignments from fast swap, before saving");
        // 	System.out.println("Parent:");
        // 	for(p = parent.last; p != null; p = p.prev){
        // 	    System.out.print(p+"("+(p.selected ? "y" : "n")+","+(p.emptyWindow ? "y" : "n")+") ");
        // 	}
        // 	System.out.println();
        // 	System.out.println("Parent.left:");
        // 	for(p = parent.last; p != null; p = p.prev){
        // 	    if(p.left != null){
        // 		System.out.print(p.left+"("+(p.left.selected ? "y" : "n")+","+(p.left.emptyWindow ? "y" : "n")+") ");
        // 	    }
        // 	    else{
        // 		System.out.print(p.left+"(n,n) ");
        // 	    }
        // 	}
        // 	System.out.println();
        // 	System.out.println("Parent.right:");
        // 	for(p = parent.last; p != null; p = p.prev){
        // 	    if(p.right != null){
        // 		System.out.print(p.right+"("+(p.right.selected ? "y" : "n")+","+(p.right.emptyWindow ? "y" : "n")+") ");
        // 	    }
        // 	    else{
        // 		System.out.print(p.right+"(n,n) ");
        // 	    }
        // 	}
        // 	System.out.println();
        // 	if(parent.left != null){
        // 	    System.out.println("Left:");
        // 	    for(p=parent.left.last; p != null; p = p.prev){
        // 		System.out.print(p+"("+(p.selected ? "y" : "n")+","+(p.emptyWindow ? "y" : "n")+") ");
        // 	    }
        // 	    System.out.println();
        // 	    System.out.println("Left parents:");
        // 	    for(p=parent.left.last; p != null; p = p.prev){
        // 		System.out.print(p.parent+" ");
        // 	    }
        // 	    System.out.println();
        // 	    System.out.println("Right:");
        // 	    for(p=parent.right.last; p != null; p = p.prev){
        // 		System.out.print(p+"("+(p.selected ? "y" : "n")+","+(p.emptyWindow ? "y" : "n")+") ");
        // 	    }
        // 	    System.out.println();
        // 	    System.out.println("Right parents:");
        // 	    for(p=parent.right.last; p != null; p = p.prev){
        // 		System.out.print(p.parent+" ");
        // 	    }
        // 	    System.out.println();

        // 	    parent.checkPointers();	    
        // 	}

        // 	/////////////////////////////////////////////////

        //saving old alignments
        copySequence();
        brother().copySequence();
        uncle.copySequence();
        parent.fullWin();
        parent.saveWin();
        grandpa.fullWin();
        grandpa.saveWin();
        if (grandpa.parent != null) {
            grandpa.copyGreatgrandpaSequence();
        }

        //NNI
        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;

        //new parent sequence
        ret += uncle.parent.alignAllWindows();
        if (Utils.DEBUG) {
        	uncle.parent.checkPointers();
        }

        // 	System.out.println("Aligned uncle.parent, pointers:");
        // 	System.out.println("Left:");
        // 	uncle.parent.left.printPointers();
        // 	System.out.println("Right:");
        // 	uncle.parent.right.printPointers();

        // 	System.out.println("RET after alignment: "+ret);

        //new grandpa
        ret += grandpa.alignAllWindows();

        if (Utils.DEBUG) {
        	grandpa.checkPointers();
        }
        // 	System.out.println("aligned grandpa, RET: "+ret);

        // 	System.out.println("Aligned grandpa, pointers:");
        // 	System.out.println("Left:");
        // 	grandpa.left.printPointers();
        // 	System.out.println("Right:");
        // 	grandpa.right.printPointers();


        //if there is a greatgrandpa, then we align it...
        if (grandpa.parent != null) {
            //System.out.println("Greatgrandpa in fastSwap: "+grandpa.parent);
            AlignColumn gg = grandpa.parent.last, ggf = gg;
            g = grandpa.last;
            gf = g;
            while (g != null) {
                gf = gf.findWindowStart();
                ggf = ggf.findWindowStart();
                grandpa.winLast = g;
                grandpa.winFirst = gf;
                grandpa.parent.winLast = gg;
                grandpa.parent.winFirst = ggf;
                if (gg.emptyWindow || true) {
                    ret += grandpa.hmm2Align();
                }
                g = gf.prev;
                gf = g;
                grandpa.parent.first = ggf;
                gg = ggf.prev;
                ggf = gg;
                if (g != null) {
                    g.parent = gg;
                    g.orphan = false;
                    if (grandpa.parent.left == grandpa) {
                        gg.left = g;
                    } else {
                        gg.right = g;
                    }
                }
            }

            //	    System.out.println("RET after aligning everything: "+ret);
            //System.out.println("Calculating Felsenstein, called from fastSwap");
            grandpa.parent.calcFelsen();
            if (Utils.DEBUG) {
	            grandpa.parent.checkPointers();
            }
            //System.out.println("Aligned greatgrandpa, pointers:");
            //System.out.println("Left:");
            //grandpa.parent.left.printPointers();
            //System.out.println("Right:");
            //grandpa.parent.right.printPointers();
            //System.out.println("and the pointers of the parent sequence...");
            //for(p = grandpa.parent.first; p != grandpa.parent.last; p = p.next){
            //		System.out.println("left: "+p.left+" right: "+p.right);
            //}


        }


        // 	ret += uncle.parent.hmm3AlignWithSave();
        // 	ret += grandpa.hmm3AlignWithSave();
        // 	//ret += uncle.parent.hmm3Align();
        // 	//ret += grandpa.hmm3Align();
        // 	if(grandpa.parent != null) {
        // 	    ret += grandpa.hmm2AlignWithSave();
        // 	    //ret += grandpa.hmm2Align();
        // 	    grandpa.parent.calcAllUp();
        // 	} else {
        // 	    grandpa.calcOrphan();
        // 	}
        calcOrphan();
        uncle.brother().calcOrphan();
        uncle.calcOrphan();
        uncle.calcAllUp();


        //And finally: the window selecting backproposals...
        //	System.out.println((grandpa.parent == starter ? "grandpa.parent == starter" : "grandpa.parent != starter"));
        ret += starter.backproposeAnchors();

        //	System.out.println("RET after backproposing anchors (final value) "+ret);

        return ret;
    }

    /** This function iteratively aligns the windows for topology change */
    double alignAllWindows() {
        double ret = 0.0;

        length = 0;
        AlignColumn p = last;// pf = last;
        AlignColumn l = left.last, lf = left.last;
        AlignColumn r = right.last, rf = right.last;
        //	System.out.println("Aligning sequences here: "+print());
        //	owner.printAllPointers();
        //System.out.println("Printing the pointers +alignments from alignAll");
        //System.out.println("Parent:");
        //for(p = this.last; p != null; p = p.prev){
        //  System.out.print(p+"("+(p.selected ? "y" : "n")+","+(p.emptyWindow ? "y" : "n")+") ");
        //}
        //System.out.println();
        //System.out.println("Parent.left:");
        //for(p = this.last; p != null; p = p.prev){
        //  if(p.left != null){
        //	System.out.print(p.left+"("+(p.left.selected ? "y" : "n")+","+(p.left.emptyWindow ? "y" : "n")+") ");
        //  }
        //  else{
        //	System.out.print(p.left+"(n,n) ");
        //  }
        //}
        //System.out.println();
        //System.out.println("Parent.right:");
        //for(p = this.last; p != null; p = p.prev){
        //  if(p.right != null){
        //	System.out.print(p.right+"("+(p.right.selected ? "y" : "n")+","+(p.right.emptyWindow ? "y" : "n")+") ");
        //  }
        //  else{
        //	System.out.print(p.right+"(n,n) ");
        //  }
        //}
        //System.out.println();
        if (left != null) {
            //    System.out.println("Left:");
            //for(p=this.left.last; p != null; p = p.prev){
            //	System.out.print(p+"("+(p.selected ? "y" : "n")+","+(p.emptyWindow ? "y" : "n")+") ");
            //}
            //System.out.println();
            //  System.out.println("Left parents:");
            //for(p=this.left.last; p != null; p = p.prev){
            //		System.out.print(p.parent+" ");
            //}
            //System.out.println();
            //System.out.println("Right:");
            //for(p=this.right.last; p != null; p = p.prev){
            //	System.out.print(p+"("+(p.selected ? "y" : "n")+","+(p.emptyWindow ? "y" : "n")+") ");
            //}
            //System.out.println();
            //System.out.println("Right parents:");
            //for(p=this.right.last; p != null; p = p.prev){
            //	System.out.print(p.parent+" ");
            //}
            // System.out.println();

        	if (Utils.DEBUG) {
        		checkPointers();
        	}
        }

        p = last; //pf = last;
        l = left.last;
        lf = left.last;
        r = right.last;
        rf = right.last;


        while (l != null) {
            //    System.out.println("****Start finding the left window: ");
            lf = lf.findWindowStart();
            //	    System.out.println("****Finished finding the left window: ");
            rf = rf.findWindowStart();
            winLength = 0;
            left.winFirst = lf;
            left.winLast = l;
            right.winFirst = rf;
            right.winLast = r;
            winLast = p;
            winFirst = p;
            if (p.emptyWindow || true) {
                //	System.out.println("Aligning by hmm3Align()!!! left first: "+left.winFirst+" left last: "+left.winLast+
                //		   " right first: "+right.winFirst+" right last "+right.winLast);
                //System.out.println("parent first: "+winFirst+" its children, left, right "+winFirst.left+", "+winFirst.right+
                //		   "parent last:  "+winLast+" its children, left, right "+winLast.left+", "+winLast.right);

                ret += hmm3Align();
                p = winFirst;
            }
            //	    else{
            //	pf.next = p;
            //p.prev = pf;
            //p = pf;
            //	System.out.println("Not aligning this!!! left first: "+left.winFirst+" left last: "+left.winLast+
            //		   " right first: "+right.winFirst+" right last "+right.winLast);
            //System.out.println("parent first: "+winFirst+" its children, left, right "+winFirst.left+", "+winFirst.right+
            //		   "parent last:  "+winLast+" its children, left, right "+winLast.left+", "+winLast.right);

            //	    }
            l = lf.prev;
            lf = lf.prev;
            r = rf.prev;
            rf = rf.prev;
            if (l != null) {
                p = new AlignColumn(p, true);
                length++;
                p.left = l;
                p.right = r;
                l.parent = p;
                l.orphan = false;
                r.parent = p;
                r.orphan = false;
                p.selected = l.selected;
                p.emptyWindow = l.emptyWindow;
            }
            first = p;
        }
        calcFelsen();

        //System.out.println("\n\n******* Printing the pointers +alignments from alignAll after aligning ****");
        //System.out.println("Parent:");
        //for(p = this.last; p != null; p = p.prev){
        //  System.out.print(p+" ");
        //}
        //System.out.println();
        //System.out.println("Parent.left:");
        //for(p = this.last; p != null; p = p.prev){
        //  System.out.print(p.left+" ");
        //}
        //System.out.println();
        //System.out.println("Parent.right:");
        //for(p = this.last; p != null; p = p.prev){
        //  System.out.print(p.right+" ");
        //	}
        //System.out.println();
        if (left != null) {
            //    System.out.println("Left:");
            //for(p=this.left.last; p != null; p = p.prev){
            //	System.out.print(p+" ");
            //}
            //System.out.println();
            //System.out.println("Left parents:");
            //for(p=this.left.last; p != null; p = p.prev){
            //	System.out.print(p.parent+" ");
            //}
            //System.out.println();
            //System.out.println("Right:");
            //for(p=this.right.last; p != null; p = p.prev){
            //	System.out.print(p+" ");
            //}
            //System.out.println();
            //System.out.println("Right parents:");
            //for(p=this.right.last; p != null; p = p.prev){
            //	System.out.print(p.parent+" ");
            //}
            //System.out.println();

        	if (Utils.DEBUG) {
        		checkPointers();
        	
	            //checking pointer integrity
	            for (AlignColumn c = left.first; c != null; c = c.next) {
	                p = first;
	                while (c.parent != p && p != null) {
	                    p = p.next;
	                }
	                if (p == null) {
	                    throw new Error("children does not have a parent!!!" + this + " " + this.print(3));
	                }
	            }
	            for (AlignColumn c = right.first; c != null; c = c.next) {
	                p = first;
	                while (c.parent != p && p != null) {
	                    p = p.next;
	                }
	                if (p == null) {
	                    throw new Error("children does not have a parent!!!" + this + " " + this.print(3));
	                }
	            }
        	}


        }

        //	left.printPointers();
        //right.printPointers();
        //System.out.println("and the pointers of the parent sequence...");
        //for(p = first; p != last; p = p.next){
        //  System.out.println("left: "+p.left+" right: "+p.right);
        //}
        //	String[] s = left.printedAlignment();
        //System.out.println("Left:\n"+s[0]+"\n"+s[1]);
        //s = right.printedAlignment();
        //System.out.println("Ritght:\n"+s[0]+"\n"+s[1]);

        //System.out.println("\n\n************ End of alignAllWindows ***************\n\n");

        return ret;
    }

    /**
     * Restores the exact state just before the call of {@link #fastSwapWithUncle()}.
     * Must be called on ex-uncle node.
     * Assumes {@code this} has a non-null grandparent.
     */
    public void fastSwapBackUncle() {
        Vertex uncle = parent.brother(), grandpa = parent.parent;

        fullWin();
        brother().fullWin();
        parent.fullWin();
        uncle.fullWin();
        grandpa.fullWin();
        if (grandpa.parent != null) {
            grandpa.parent.fullWin();
        }


        parentNewChild(uncle);                    // order is important!
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;


        //	System.out.println("Starting alignment restauration...");
        grandpa.alignRestore();
        //System.out.println("End alignment restauration...");
    }

    /**
     * This function selects anchor homologous columns that are not changed during topology change
     * returns with the -logarithm of proposal probability!!!
     */
    double selectAnchors() {
        double ret = 0.0;

        AlignColumn actual = first;
        while (actual != null) {
            if (actual.isHomologous()) {
                if (Utils.generator.nextDouble() < SELECTING) {
                    actual.selectDown();
                    ret -= Math.log(SELECTING);
                    // System.out.println("Homologous, selected!");
                    if (actual.isEmptyWindow()) {
                        if (Utils.generator.nextDouble() < EMPTY_WINDOW) {
                            actual.setEmptyWindowDown(true);
                            ret -= Math.log(EMPTY_WINDOW);
                        } else {
                            actual.setEmptyWindowDown(false);
                            ret -= Math.log(1.0 - EMPTY_WINDOW);
                        }
                    } else {// we need this since later on we will not know if it is not an empty window...
                        actual.setEmptyWindowDown(true);
                    }
                } else {
                    actual.selected = false;
                    ret -= Math.log(1.0 - SELECTING);
                    // System.out.println("Homologous, not selected!");
                }
            } else {
                actual.selected = false;
                //	System.out.println("Not homologous!");
            }
            actual = actual.next;
            //System.out.println("...");
        }

        return ret;
    }

    /**
     * This function selects anchor homologous columns that are not changed during topology change
     * returns with the logarithm of proposal probability!!!
     */
    double backproposeAnchors() {
        double ret = 0.0;

        AlignColumn actual = first;
        //	while(actual.next != null){
        while (actual != null) {
            if (actual.isHomologous()) {
                if (actual.selected) {
                    ret += Math.log(SELECTING);
                    //    System.out.println("Homologous, selected!");
                    if (actual.isEmptyWindow()) {
                        if (actual.emptyWindow) {
                            ret += Math.log(EMPTY_WINDOW);
                        } else {
                            ret += Math.log(1.0 - EMPTY_WINDOW);
                        }
                    }
                } else {
                    ret += Math.log(1.0 - SELECTING);
                    //    System.out.println("Homologous, not selected!");
                }
            }
            //    else{
            //	System.out.println("It is not homologous!!!");
            //    }
            actual = actual.next;
            //System.out.println("...");
        }


        //for(actual = first; actual != null; actual = actual.next){
        //  System.out.println("left: "+actual.left+" right: "+actual.right);
        //}

        //	String[] s = printedMultipleAlignment();
        //for(int i = 0; i < s.length; i++){
        //  System.out.println(s[i]);
        //}

        return ret;
    }

    /** this function set selected to false for all AlignColumns */
    void setAllAlignColumnsUnselected() {
        AlignColumn a = first;
        while (a != null) {
            a.selected = false;
            a = a.next;
        }
    }

    /** This function copies the sequence to the old vertex */
    void copySequence() {

        fullWin();
        saveWin();
        AlignColumn a = last;
        AlignColumn prev = last;
        while (a.prev != null) {
            a = a.prev;
            prev = new AlignColumn(prev, false);
            prev.saveBoth(a);
        }
        winFirst = prev;
        first = prev;

        if (left != null) {
            //winFirst = first;
            //winLast = last;
            //old.winFirst = old.first;
            //old.winLast = old.last;
            left.toggleUp(true);
            //System.out.println("Left is ready");
            right.toggleUp(true);
        }
    }

    /**
     * This function is used for special copying the geatgrandpa's sequence in topology change.
     * The function is called to the grandpa, whose parent is the greatgrandpa
     */
    void copyGreatgrandpaSequence() {
        parent.fullWin();
        parent.saveWin();
        AlignColumn a = parent.last;
        AlignColumn prev = parent.last;
        boolean isLeft = (parent.left == this);
        while (a.prev != null) {
            a = a.prev;
            prev = new AlignColumn(prev, true);
            if (isLeft) {
                prev.saveRight(a);
            } else {
                prev.saveLeft(a);
            }
            prev.saveParent(a);
        }
        parent.winFirst = prev;
        parent.first = prev;
        if (parent.left == this) {
            parent.right.toggleUp(true);
            //	    System.out.println("Toggling up right!!!");
        } else {
            parent.left.toggleUp(true);
            //System.out.println("Toggling up left!!!");
        }
        //	System.out.println("Finished copying the greatgrandpa, pointers:\nLeft:");
        //parent.left.printPointers();
        //System.out.println("Right:");
        //parent.right.printPointers();
        //System.out.println("greatgrandpa("+parent+") pointers: ");
        //for(a = parent.first; a != parent.last; a = a.next){
        //  System.out.println("Left: "+a.left+" right: "+a.right);
        //}


    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /** Calculates Felsenstein and indel likelihoods up to root, starting from `parent' */
    public void calcAllUp() {
        for (Vertex v = parent; v != null; v = v.parent) {
            v.calcFelsen();
            v.calcOrphan();
            v.calcIndelLogLike(false);
        }
    }

    void parentNewChild(Vertex child) {
        child.last.parent = parent.last;
        if (parent.left == this) {
            parent.left = child;
            parent.last.left = child.last;
        } else {
            parent.right = child;
            parent.last.right = child.last;
        }
    }

    /**
     * Returns the alignment of this {@link Vertex} to its parent as an integer array:
     * for each character position <i>i</i> in this sequence it tells the character index
     * in the parent that <i>i</i> is aligned to or <i>-j-1</i> if it is not aligned
     * (gap: -) and the next character in the parent has index <i>j</i> (indices start from 0).
     * If this vertex does not have a parent an all-zero array is returned that has the same
     * length as the sequence in this vertex.
     */
    int[] getAlign() {
        int[] align = new int[length];
        Vertex vp = parent;
        if (vp != null) {
            AlignColumn c = first, p = vp.first;
            int cn = 0, pn = 0;
            while (c != last || p != vp.last) {
                if (c.parent != p) {            // deletion (* -)
                    pn++;
                    p = p.next;
                } else if (c.orphan) {        // insertion (- *)
                    align[cn] = -pn - 1;
                    cn++;
                    c = c.next;
                } else {                    // substitution (* *)
                    align[cn] = pn;
                    pn++;
                    p = p.next;
                    cn++;
                    c = c.next;
                }
            }
        }
        return align;
    }

    /** Returns the Felsenstein loglikelihood vector for the sequence in this vertex. */
    double[][] getFelsen() {
        double[][] felsen = new double[length][];

        int i = 0;
        for (AlignColumn c = first; c != last; c = c.next) {
            felsen[i++] = c.seq.clone();
        }
        return felsen;
    }
    
    private boolean inWindow = false;

    String[] printedAlignmentInWindow() {
    	inWindow = true;
    	String[] s = printedAlignment();
    	inWindow = false;
    	return s;
    }
    public String[] printedAlignment() {
        //	try{
        //	printPointers();
        if (parent == null) {
            String[] s = new String[2];
            AlignColumn a = inWindow ? winFirst : first;
            //AlignColumn a = winFirst;
            s[0] = "";
            s[1] = "";
            while (a != (inWindow ? winLast : last)) {
            //while (a != winLast) {
                s[0] += "-";
                s[1] += a.mostLikely();
                a = a.next;
            }
            return s;
        }
        //	System.out.println("The parent !=null!!!");
        String[] s = new String[2];
        s[0] = "";
        s[1] = "";
        AlignColumn p = inWindow ? parent.winFirst : parent.first;
        //AlignColumn p = parent.winFirst;
        AlignColumn c = inWindow ? winFirst : first;
        //AlignColumn c = winFirst;
        while (p != (inWindow ? parent.winLast : parent.last) || c != (inWindow ? winLast : last)) {
        //while (p != parent.winLast || c != winLast) {
            // System.out.println(".");
            if (c.parent != p) {
                if (c.prev == null || (c.prev.parent != p || c.prev.orphan)) {
                    //	s[0] += "*";
                    s[0] += p.mostLikely();
                    s[1] += "-";
                }
                p = p.next;
                //	System.out.println("We increased p");
            } else {
                if (c.orphan) {
                    s[0] += "-";
                } else {
                    //	    s[0] += "*";
                    //	double max = 0.0;
                    //	char character = '*';
                    //	for(int i = 0; i < p.seq.length; i++){
                    //		if(p.seq[i] > max){
                    //			character =  owner.substitutionModel.alphabet[i];
                    //			max = p.seq[i];
                    //			//			s[1] += owner.substitutionModel.alphabet[i];
                    //			//found = true;
                    //		}
                    //	}
                    s[0] += p.mostLikely();
                }
                //	boolean found = false;
                //double max = 0.0;
                //char character = '*';
                //for(int i = 0; !found && i < c.seq.length; i++){
                //	// if(c.seq[i] > 0.9999){
                //	if(c.seq[i] > max){
                //		character =  owner.substitutionModel.alphabet[i];
                //		max = c.seq[i];
                //		//			s[1] += owner.substitutionModel.alphabet[i];
                //		//found = true;
                //	}
                //}
                //if(!found){
                //  s[1] += "*";
                //	}
                s[1] += c.mostLikely();
                c = c.next;
                //	System.out.println("We increased c");
            }
        }

        //System.out.println("Alignment here: "+print()+"\n"+s[0]+"\n"+s[1]);

        return s;
        //	}
        /*	
		catch(NullPointerException e){
		System.out.println(print());
		AlignColumn c = first;
		AlignColumn p = parent.first;
		while(p != null){
		System.out.print(p+" ");
		p = p.next;
		}
		System.out.println();
		while(c != null){
		System.out.print(c.parent+" ");
		c = c.next;
		}
		System.out.println();
		throw new Error("Couldn't find the parent of one of the nodes!");

		}
		 */
        //	return new String[] {" "," "};
    }

    //////////////////////////////////////////////////
    String[] printedMultipleAlignment() {
        if (left == null && right == null) {
            String n = new String(name);
            /*
			while(n.length() < 30){
				n += " ";
			}
			if(n.length() > 30){
				n = (String)n.subSequence(0,30);
			}
			 */
            String s = "";
            String s1 = printedAlignment()[1];
            for (int i = 0; i < s1.length(); i++) {
                if (s1.charAt(i) != '-') {
                    s += s1.charAt(i);
                }
            }
            return new String[]{n + "\t" + s};
        }
        String[] d1 = left.printedMultipleAlignment();
        String[] d2 = right.printedMultipleAlignment();

        String[] e1 = left.printedAlignment();
        String[] e2 = right.printedAlignment();

        //		System.out.println("We are in subtree: "+print());
        //		System.out.println("d1:");
        //		for(int i = 0; i < d1.length; i++){
        //		System.out.println(d1[i]);
        //		}
        //		System.out.println("d2:");
        //		for(int i = 0; i < d2.length; i++){
        //		System.out.println(d2[i]);
        //		}
        //		System.out.println("e1:");
        //		for(int i = 0; i < e1.length; i++){
        //		System.out.println(e1[i]);
        //		}
        //		System.out.println("e2:");
        //		for(int i = 0; i < e2.length; i++){
        //		System.out.println(e2[i]);
        //		}
        //		System.out.println("pointers:");
        //		left.printPointers();
        //		right.printPointers();

        //		System.out.println("Starting the calculation...");
        int tabPosition = d1[0].indexOf('\t');
        String[] s = new String[d1.length + d2.length + 1];
        s[0] = "";
        for (int i = 0; i < tabPosition; i++) {
            s[0] += " ";
        }
        s[0] += "\t";
        for (int i = 0; i < d1.length; i++) {
            s[i + 1] = (String) d1[i].subSequence(0, tabPosition + 1);
        }
        for (int i = d1.length + 1; i < s.length; i++) {
            s[i] = (String) d2[i - d1.length - 1].subSequence(0, tabPosition + 1);
        }
        int x = 0;
        int x1 = tabPosition + 1;
        int y = 0;
        int y1 = tabPosition + 1;
        //	System.out.println("before the while cycle...");
        while (x < e1[0].length() || y < e2[0].length()) {
            while (x < e1[0].length() && e1[0].charAt(x) == '-') {
                while (d1[0].charAt(x1) == '-') {
                    s[0] += "-";
                    for (int i = 0; i < d1.length; i++) {
                        s[i + 1] += d1[i].charAt(x1);
                    }
                    for (int i = d1.length + 1; i < s.length; i++) {
                        s[i] += "-";
                    }
                    x1++;
                }
                s[0] += "-";
                for (int i = 0; i < d1.length; i++) {
                    s[i + 1] += d1[i].charAt(x1);
                }
                for (int i = d1.length + 1; i < s.length; i++) {
                    s[i] += "-";
                }
                x++;
                x1++;
            }
            while (y < e2[0].length() && e2[0].charAt(y) == '-') {
                while (d2[0].charAt(y1) == '-') {
                    s[0] += "-";
                    for (int i = 1; i < d1.length + 1; i++) {
                        s[i] += "-";
                    }
                    for (int i = d1.length + 1; i < s.length; i++) {
                        s[i] += d2[i - d1.length - 1].charAt(y1);
                    }
                    y1++;
                }
                s[0] += "-";
                for (int i = 1; i < d1.length + 1; i++) {
                    s[i] += "-";
                }
                for (int i = d1.length + 1; i < s.length; i++) {
                    s[i] += d2[i - d1.length - 1].charAt(y1);
                }
                y++;
                y1++;
            }
            if (x < e1[0].length() && y < e2[0].length()) {
                while (x1 < d1[0].length() && d1[0].charAt(x1) == '-') {
                    s[0] += "-";
                    for (int i = 0; i < d1.length; i++) {
                        s[i + 1] += d1[i].charAt(x1);
                    }
                    for (int i = d1.length + 1; i < s.length; i++) {
                        s[i] += "-";
                    }
                    x1++;
                }
                while (y1 < d2[0].length() && d2[0].charAt(y1) == '-') {
                    s[0] += "-";
                    for (int i = 1; i < d1.length + 1; i++) {
                        s[i] += "-";
                    }
                    for (int i = d1.length + 1; i < s.length; i++) {
                        s[i] += d2[i - d1.length - 1].charAt(y1);
                    }
                    y1++;
                }
                //	s[0] += "*";
                s[0] += e1[0].charAt(x);
                for (int i = 0; i < d1.length; i++) {
                    s[i + 1] += (e1[1].charAt(x) != '-' ? d1[i].charAt(x1) : "-");
                }
                x1 += (e1[1].charAt(x) != '-' ? 1 : 0);
                x++;
                for (int i = d1.length + 1; i < s.length; i++) {
                    s[i] += (e2[1].charAt(y) != '-' ? d2[i - d1.length - 1].charAt(y1) : "-");
                }
                y1 += (e2[1].charAt(y) != '-' ? 1 : 0);
                y++;
            }
        }
        for (x1 = x1 + 0; x1 < d1[0].length(); x1++) {
            s[0] += "-";
            for (int i = 0; i < d1.length; i++) {
                s[i + 1] += d1[i].charAt(x1);
            }
            for (int i = d1.length + 1; i < s.length; i++) {
                s[i] += "-";
            }
        }
        for (y1 = y1 + 0; y1 < d2[0].length(); y1++) {
            s[0] += "-";
            for (int i = 0; i < d1.length; i++) {
                s[i + 1] += "-";
            }
            for (int i = d1.length + 1; i < s.length; i++) {
                s[i] += d2[i - d1.length - 1].charAt(y1);
            }
        }

        //		System.out.println("End of while cycle");


        //		System.out.println("The alignment contains "+s.length+" column");
        //		for(int i = 0; i < s.length; i++){
        //		System.out.println(s[i]);
        //		}
        //		System.out.println();

        return s;
    }

    void printPointers() {
        AlignColumn c = first;
        AlignColumn p = parent.first;
        while (p != null) {
            System.out.print(p + " ");
            p = p.next;
        }
        System.out.println();
        while (c != null) {
            System.out.print(c.parent + " ");
            c = c.next;
        }
        System.out.println("\n");

        c = first;
        p = parent.first;
        while (p.next != null) {
            System.out.print(p.mostLikely() + "");
            p = p.next;
        }
        System.out.println();
        while (c != null) {
            if (c.parent.next != null) {
                System.out.print(c.parent.mostLikely() + "");
            }
            c = c.next;
        }
        System.out.println();
        c = first;
        while (c != null) {
            if (c.parent.next != null) {
                System.out.print((c.orphan ? "y" : "n"));
            }
            c = c.next;
        }
        System.out.println("\n");


    }
    
    void collectIndicesBelow(Vertex v, ArrayList<Integer> indices) {
    	if(v == null) return;
    	indices.add(v.index);
    	collectIndicesBelow(v.left, indices);
    	collectIndicesBelow(v.right, indices);
	}

    public void printToScreenAlignment(int windowStart, int windowEnd) {
		int nGapsInWindow = 0;
    	if (parent != null) {
    		String[] s0 = printedAlignmentInWindow();
    		System.out.println("parent\t"+s0[0]);
    		System.out.println("child\t"+s0[1]+"\n");
//    		for (int k=0; k< windowEnd; k++) {
//    			nGapsInWindow += (s0[1].charAt(k) == '-') ? 1 : 0;
//    		}
    	}
		String[] s = owner.getState().getFullAlign();
		int start = 0, end = 0;
		int i = 0, j = 0;
		// The stuff below is not correct -- window may include gaps on `this'
		while (i < windowStart) {
			while (s[index].charAt(start++) == '-') { }
			++i;
		}
		while (j < windowEnd) {
			while (s[index].charAt(end++) == '-') { }
			++j;
		}
		ArrayList<Integer> indices = new ArrayList<Integer>(); 
//		collectIndicesBelow(this,indices);
//    	for (int ind : indices) {
//    		System.out.println(owner.vertex[ind].name+"\t"+s[ind].substring(start,end));
//    	}
    	for (int ii=0; ii<s.length; ii++) {
    		System.out.println(owner.vertex[ii].name+"\t"+s[ii]);
    	}
    }
    

    /** Returns most likely sequence of this vertex or the original sequence for leaves. */
    String sequence() {
    	if(seq != null)
    		return seq;
        StringBuilder s = new StringBuilder();
        for (AlignColumn a = first; a != last; a = a.next)
            s.append(a.mostLikely());
        return s.toString();
    }


    /**
     * Calculate the sum of branch lengths below on this vertex.
     * Used to generate the prior for the whole state.
     * @return sum of branchlengths below this node.
     */
    public double calcSumOfEdges() {
        if (left == null) {
            return edgeLength;
        } else return edgeLength + left.calcSumOfEdges() + right.calcSumOfEdges();
    }


    /**
     * Calculate the maximum depth below on this vertex.
     * Used in tree visualisation.
     * @return maximum depth below this vertex.
     */
    public double maxDepth() {
        if (left == null) {
            return edgeLength;
        } else return edgeLength + Math.max(left.maxDepth(), right.maxDepth());
    }

    void compareLike(Vertex vert) {
        if (left != null && right != null && vert.left != null && vert.right != null) {
            left.compareLike(vert.left);
            right.compareLike(vert.right);
        }
        System.out.println("this.indel=" + indelLogLike + " tree.indel=" + vert.indelLogLike);
        System.out.println("indel " + (Math.abs(indelLogLike - vert.indelLogLike) < 1e-8 ? "ok" : "error"));
        System.out.println("this.orphan=" + orphanLogLike + " tree.orphan=" + vert.orphanLogLike);
        System.out.println("orphan " + (Math.abs(orphanLogLike - vert.orphanLogLike) < 1e-8 ? "ok" : "error"));
    }


    /** this function checks if the pointers are all right... */
    public void checkPointers() {
        //parent
        for (AlignColumn p = first; p != null; p = p.next) {
            if (p.left != null && (p.left.orphan || p.left.parent != p)) {
                throw new Error("Problem is vertex " + this + ":\np is: " + p + " p.left is " + p.left + " p.left.orphan: " + p.left.orphan +
                        " p.left.parent: " + p.left.parent);
            }
            if (p.right != null && (p.right.orphan || p.right.parent != p)) {
                throw new Error("Problem is vertex " + this + ":\np is: " + p + " p.right is " + p.right + " p.right.orphan: " + p.right.orphan +
                        " p.right.parent: " + p.right.parent);
            }
        }
        for (AlignColumn l = left.first; l != null; l = l.next) {
            if (!l.orphan) {
                if (l.parent == null || l.parent.left != l) {
                    throw new Error("Problem in vertex " + this + ":\nl is: " + l + (l.parent == null ? " l does not have a parent" : " l parent is: " + l.parent + " l parent left is: " + l.parent.left));
                }
            }
        }
        for (AlignColumn r = right.first; r != null; r = r.next) {
            if (!r.orphan) {
                if (r.parent == null || r.parent.right != r) {
                    throw new Error("Problem in vertex " + this + ":\nr is: " + r + (r.parent == null ? " r does not have a parent" : " r parent is: " + r.parent + " r parent right is: " + r.parent.right));
                }
            }
        }
    }

//    /**
//     * This function is merely for testing/debugging purposes
//     * @param args Arguments are not used, all input is directly written into the
//     *             function.
//     */
//    public static void main(String[] args) throws IOException {
//        try {
//            Tree tree = new Tree(new String[]{"qqqqqqqqqqqqkkkwwwwwlidwwwwwkkk",
//                    "kkkwwwwwlidwwwwwkkk",
//                    "qqqqqqqqqqqqkkkwwwwwlidwwwwwkkk",
//                    "kkkwwwwwlidwwwwwkkkeqeq",
//                    "kkkwwwwwlidwwwwwkkkddkldkl",
//                    "kkkwwwwwlidwwwwwkkkeqiqii",
//                    "kkkwwwwwlidwwwwwkkkddkidkil",
//                    "kkkwwwwwlidwwwwwkkkeqiq",
//                    "kkkwwwwwlidwwwwwkkkddkldkll",
//                    "kkkwwwwwlidwwwwwkkkddkldkil"},
//                    new String[]{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"},
//                    new Dayhoff(),
//                    new Blosum62(), "");
//            System.out.println(tree.printedTree());
//            tree.root.calcFelsRecursively();
//            System.out.println(tree.root.orphanLogLike);
//
////			String[] s = tree.printedAlignment("StatAlign");
//            /*String[] s = null;
//            for (int i = 0; i < s.length; i++) {
//                System.out.println(s[i]);
//            }*/
//
//            /////////////////////////////////
//            double[] weights = new double[tree.vertex.length];
//            tree.countLeaves(); // calculates recursively how many leaves we have below this node
//            for (int j = 0; j < 100; j++) {
//                for (int i = 0; i < weights.length; i++) {
//                    weights[i] = Math.pow(tree.vertex[i].leafCount, Mcmc.LEAFCOUNT_POWER);
//                }
//                int k = Utils.weightedChoose(weights, null);
//                tree.vertex[k].selectSubtree(Mcmc.SELTRLEVPROB, 0);
//                System.out.println("Selected vertices: ");
//                for (int i = 0; i < tree.vertex.length; i++) {
//                    if (tree.vertex[i].selected) {
//                        System.out.print(i + " ");
//                        tree.vertex[i].selected = false;
//                    }
//                }
//                System.out.println("\n");
//
//            }
//        } catch (StoppedException e) {
//            // stopped during tree construction
//        }
//
//    }

}
