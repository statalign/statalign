package statalign.base;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.FieldPosition;
import java.util.ArrayList;
import static statalign.base.AlignColumn.*;

import statalign.ui.ErrorMessage;

/**
 * This is a vertex of the tree.
 * The 'hardcore' functions are implemented in this class, developers are suggested
 * not change functions in it. The implemented functions are quite unreadable, since we
 * opted for efficiency and not readability. You should be able to develop novel
 * functionality of the software package (postprocessing, substitution models, etc.)
 * without touching this class.
 * @author miklos, novak, herman
 */
public class Vertex {

	/** The number of digits in rounding when printing the tree. */
    static final int ROUNDING = 100; 
    /** The probability for selecting a homologous column not to be 
     *  changed during topology changing 
     */
    static final double SELECTING = 0.5; 
    static int DEBUG = 0;
    static final double PROPOSAL_HEAT = 1.0;
    /** The probability that an empty window will be realigned*/
    static final double EMPTY_WINDOW = 0.01; 
    
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
        //double heat = PROPOSAL_HEAT;
        double heat = 1.0;
        hmm2PropTransMatrix = new double[hmm2TransMatrix.length][];
        for (int i = 0; i < hmm2TransMatrix.length; i++) {
            hmm2PropTransMatrix[i] = hmm2TransMatrix[i].clone();
            /**/
			double tempSum = Utils.log0;
			for (int j = 0; j < hmm2PropTransMatrix[i].length; j++) {

				hmm2PropTransMatrix[i][j] = hmm2PropTransMatrix[i][j]* heat;
				//hmm2PropTransMatrix[i][j] = -Math.pow(Math.abs(hmm2PropTransMatrix[i][j]),heat);

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

            
			//double heat = PROPOSAL_HEAT;
            double heat = 1.0;
			
//			double heat = owner.heat;
//			
			for (int i = 0; i < hmm3RedTransMatrix.length; i++) {
				double tempSum = Utils.log0;
				for (int j = 0; j < hmm3RedTransMatrix[i].length; j++) {
					hmm3RedTransMatrix[i][j] = hmm3RedTransMatrix[i][j] * heat;
					//hmm3RedTransMatrix[i][j] = -Math.pow(Math.abs(hmm3RedTransMatrix[i][j]),heat);
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
        	b.append(name.replaceAll(" ", "")+"["+index+"]");
        	b.append(':');
        	printDf.format(edgeLength, b, new FieldPosition(1));
        } else {
        	b.append('(');
        	left.print(b);
        	b.append(',');
        	right.print(b);
        	b.append(')'+"["+index+"]");
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

    public void recomputeCheckLogLike() {
        if (left != null && right != null) {
            // System.out.println("calling the left child");
            left.recomputeCheckLogLike();
            //System.out.println("calling the right child");
            right.recomputeCheckLogLike();
            //calcUpperWithCheck(); 
            // we shouldn't require all upp vectors to be
            // recomputed after every move, because that's unnecessary
        }
        calcFelsenWithCheck();
        calcOrphanWithCheck();
        calcIndelLogLikeWithCheck();
    }

    /**
     * This function calculates the upper probability vectors, which contain
     * the partial likelihoods for everything except this subtree.
     * Requires the felsenstein partial likelihoods to have already been 
     * computed for everything below.
     */
    public void calcUpperRecursively() {
    	calcUpper();
    	//updateTransitionMatrix();
    	if (left != null && right != null) {
            // System.out.println("calling the left child");
            left.calcUpperRecursively();
            //System.out.println("calling the right child");
            right.calcUpperRecursively();
        }
    }
    public void calcUpperRecursivelyWithCheck() {
    	calcUpperWithCheck();
    	if (left != null && right != null) {
            // System.out.println("calling the left child");
            left.calcUpperRecursivelyWithCheck();
            //System.out.println("calling the right child");
            right.calcUpperRecursivelyWithCheck();
        }
    }
    void calcUpper(boolean withCheck, boolean withoutBrother) {
    	AlignColumn v; 
    	Vertex brother = null;
    	    	
    	boolean isLeft = false;
    	
    	if (this != owner.root) {
    		brother = brother();
    		isLeft = (this == parent.left);
    		if (Utils.DEBUG) {
    			parent.calcUpperWithCheck();
    			brother.calcFelsenWithCheck();
    		}
    	}   	
    	
        for (v = first; v != last; v = v.next) {
        	double[] upp_old = null;
        	if (withCheck) {
        		upp_old = new double[v.upp.length];
        		for (int i=0; i<upp_old.length; i++) {
        			upp_old[i] = v.upp[i];
        		}
        	}
        	if (this == owner.root || v.orphan) {
        		v.upp = owner.substitutionModel.e.clone();
        	}
        	else {
        		AlignColumn p = v.parent;
        		AlignColumn b = isLeft ? p.right : p.left;
        		double[] felsen = new double[ owner.substitutionModel.e.length];       		
        		if (!withoutBrother) {
                    Utils.calcFelsen(felsen, null, null, b != null ? b.seq : null, brother.charTransMatrix);
        		}
        		else for (int i=0; i<felsen.length; i++)  felsen[i] = 1.0;
        		
        		v.upp = new double[ owner.substitutionModel.e.length];
        		for (int i=0; i<v.upp.length; i++) {
        			for (int j=0; j<v.upp.length; j++) {
        				v.upp[j] += felsen[i] * p.upp[i] * charTransMatrix[i][j];
	        		}
        		}
        	}
        	if (withCheck) {
        		boolean match = true;
        		for(int i = 0; i < upp_old.length && match; i++){
        			match = Math.abs(upp_old[i]/v.upp[i] - 1.0) < 0.0001;
        		}
        		if (!match) {
        			String info = "";
        			for(int i = 0; i < upp_old.length; i++){
        				info += (upp_old[i])+"/"+(v.upp[i])+" ";
        			}
        			info += "\n";
        			throw new RuntimeException("Upper probabilities do not match:\n"+info);
        		}
        	}
        }
    }
    void calcUpper() {
    	calcUpper(false,false);
    }
    void calcUpperWithCheck() {
    	calcUpper(true,false);
    }
    void calcUpperWithoutBrother() {
    	calcUpper(false,true);
    }
        

    /**
     * This function calculates the Felsenstein likelihood for the subtree below <code>this</code>. 
     * When called from the root, the result is stored in <code>orphanLogLike</code> at the root.
     */
    public void calcFelsenRecursively() {
        if (left != null && right != null) {
            // System.out.println("calling the left child");
            left.calcFelsenRecursively();
            //System.out.println("calling the right child");
            right.calcFelsenRecursively();
        }
        calcFelsen();
        calcOrphan();
    }
    /** 
     * Calculates Felsenstein likelihoods of `this'. 
     * @param withCheck <code>true</code> if we are just running this function to 
     * check whether the current stored partial likelihood vectors are consistent.
     * @param ignoreChild If 0, then both children are considered (if they exist). If
     * 1, then the left child is ignored, and if 2 then the right child is ignored. 
     * This is for use during topology switching moves.
     */
    void calcFelsen(boolean withCheck, int ignoreChild) {
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
                if (withCheck) {
                	 boolean match = Utils.calcFelsenWithCheck(p.seq, fel1, left.charTransMatrix, fel2, right.charTransMatrix);
                     if (!match) {
                         throw new RuntimeException("Felsenstein does not match!");
                     }
                }
                else { 
                	if (ignoreChild == 0) {
                		Utils.calcFelsen(p.seq, fel1, left.charTransMatrix, fel2, right.charTransMatrix);
                	}
                	else if (ignoreChild == 1) { // ignore left child
                		Utils.calcFelsen(p.seq, null, null, fel2, right.charTransMatrix);
                	}
                	else if (ignoreChild == 2) { // ignore right child
                		Utils.calcFelsen(p.seq, fel1, left.charTransMatrix, null, null);
                	}
                	else throw new IllegalArgumentException("Argument `ignoreChild' must take value 0, 1 or 2.");                			
            	}
            }
        }
    }
    void calcFelsen() {
    	calcFelsen(false,0);
    }
    void calcFelsenWithCheck() {
    	calcFelsen(true,0);
    }
    /**
     * Computes the Felsenstein likelihood for this, ignoring everything below the 
     * specified child Vertex (i.e. only considering one side of the descendant tree).
     * @param child One of the two children of the current Vertex
     */
    void calcFelsenWithout(Vertex child) {
    	if (child == left) {
    		calcFelsen(false,1);
    	}
    	else if (child == right) {
    		calcFelsen(false,2);
    	}
    	else throw new IllegalArgumentException("Vertex "+child.index+" is not a child of "+index+".");
    }
    
    /**
     * Calculates the sum of orphan likelihoods in the (inclusive) subtree of `this'
     * which is then stored in `this.orphanLogLikge' (this implies that alignment to
     * parent must exists at the time of calling)
     * Saves previous orphan likelihood into `old' so it shouldn't be called twice in a row.
     * (But does not save previous Felsensteins of `this' in AlignColumn's `seq'.)
     */
    void calcOrphan(boolean withCheck) {
        old.orphanLogLike = orphanLogLike;

        //orphan likelihoods
        orphanLogLike = 0.0;

        for (AlignColumn actual = first; actual != last; actual = actual.next) {
            if (actual.parent == null || actual.orphan) {
            	orphanLogLike += Math.log(Utils.calcEmProb(actual.seq, owner.substitutionModel.e));
            }
        }

        if (left != null && right != null) {
            orphanLogLike += left.orphanLogLike + right.orphanLogLike;
        }
        if (withCheck) {
	        if (old.orphanLogLike - orphanLogLike > 1e-6) {
	            throw new RuntimeException("Problem with orphan loglike at vertex "+index+": old: " + old.orphanLogLike + " new: " + orphanLogLike);
//	            int index = 0;
//	            while (owner.vertex[index] != this) {
//	                index++;
//	            }
	        }
        }
    }
    void calcOrphan() {
    	calcOrphan(false);
    }
    void calcOrphanWithCheck() {
    	calcOrphan(true);
    }
    public void calcOrphanRecursively() {
        if (left != null && right != null) {
            left.calcOrphanRecursively();
            right.calcOrphanRecursively();
        }
        calcOrphan();
    }
    void calcOrphanRecursivelyWithCheck() {
        if (left != null && right != null) {
            left.calcOrphanRecursivelyWithCheck();
            right.calcOrphanRecursivelyWithCheck();
        }
        calcOrphanWithCheck();
    }


    public void calcIndelLogLikeRecursively() {
        if (left != null && right != null) {
            left.calcIndelLogLikeRecursively();
            right.calcIndelLogLikeRecursively();
        }
        calcIndelLogLike();
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
    void calcIndelLogLike() {
    	calcIndelLogLike(false);
    }
    void calcIndelLogLikeWithCheck() {
    	calcIndelLogLike(true);
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

        if (DEBUG==2) {
        	System.out.println("index = "+index);
        }
        //System.out.println("--------------------------------------------------");
        //printPointers();

        if (parent != null) {
            AlignColumn c = first, p = parent.first;
            int prevk = START, k;

            System.out.println(last.toString()+" "+parent.last.toString()+"\n");
            while (c != last || p != parent.last) {     
            	if (DEBUG==2)  { 
            		System.out.print(c.toString()+" "); 
            		System.out.print(c.parent.toString()+" "); 
            		System.out.println(p.toString()); 
            	}      	
                if (c.parent != p) {                            // deletion (* -), pattern code 2
                    k = emitPatt2State[2];
                    if (DEBUG==2) {
                    	System.out.print("(*-) ");
                		if (p != parent.last) System.out.println(String.valueOf(p.mostLikely())+"-");
                		else System.out.println();
                	}
                    p = p.next;                   
                } else if (c.orphan) {                        // insertion (- *), pattern code 1
                    k = emitPatt2State[1];
                    if (DEBUG==2) {
                    	System.out.print("(-*) ");
                		if (c != last) System.out.println("-"+String.valueOf(c.mostLikely()));
                		else System.out.println();
                	}
                    c = c.next;                   
                } else {                                                // substitution (* *), pattern code 3
                    k = emitPatt2State[3];
                    if (DEBUG==2) {
                    	System.out.print("(**) ");
                    	if (c != last && p != parent.last) {
                		System.out.println(String.valueOf(p.mostLikely())+String.valueOf(c.mostLikely()));
                    	}
                    	else {
                    		System.out.println(String.valueOf("-"+c.mostLikely()));
                    	}
                	}
                    p = p.next;
                    c = c.next;                   
                }
                indelLogLikeUp += hmm2TransMatrix[prevk][k];
                prevk = k;
            }

            indelLogLikeUp += hmm2TransMatrix[prevk][END];
        }
        if (DEBUG==2) System.out.println();
        
        return indelLogLikeUp;
    }

    /**
     * NB hmm3 currently does not use information about the upper parts of the tree,
     * because otherwise the recursive realignment moves would involve back-proposal
     * probabilities that would be impossible to compute. It should be possible to
     * use information about the parts of the tree above the root of the subtree,
     * which would involve repeatedly calling calcUpperWithoutBrother on the
     * downward recursion in doRecAlign, before the back-proposal probabilities are computed.
     * @return
     */
    private double[][][] hmm3ProbMatrix() {
    	return hmm3ProbMatrix(1.0);
    }
    private double[][][] hmm3ProbMatrix(double heat) {
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
                        probMatrix[i][j][k] = heat*(Math.log(emissionProb) + tr);
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
        calcIndelLogLike();
        //		System.out.println("printing the alignments at the root of this subtree: "+print());
        //		System.out.println("pointers");
        //	String[] s = left.printedAlignment();
        //System.out.println(s[0]+"\n"+s[1]+"\n");
        //s = right.printedAlignment();
        //System.out.println(s[0]+"\n"+s[1]+"\n");

        return retVal.value;
    }

    private double[][][] hmm2ProbMatrix() {
    	return hmm2ProbMatrix(1.0);
    }
    private double[][][] hmm2ProbMatrix(double heat) {
        double[] equDist = owner.substitutionModel.e;
        int hmm2Parent[] = owner.hmm2.getStateEmit()[0];
        int hmm2Child[] = owner.hmm2.getStateEmit()[1];
        int parentLen = parent.winLength, childLen = winLength;
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();
        Vertex brother = brother();

        double probMatrix[][][];                                          // DP matrix used for 2-seq HMM alignment
        probMatrix = new double[parentLen + 1][childLen + 1][END];        // don't reserve space for end state

        double emissionProb, felsen[] = new double[equDist.length], tr;
        AlignColumn c = null;                // child
        AlignColumn p = null;                // parent
        AlignColumn b;                            // brother
        int i, j, k, previ, prevj, prevk;

        /* i: parent prefix length, j: child prefix length, k: state  */
        for (i = 0; i <= parentLen; i++) {
            for (j = 0; j <= childLen; j++) {
                probMatrix[i][j][START] = (i == 0 && j == 0) ? 0.0 : Utils.log0;        // update matrix for start state
                for (k = START + 1; k < END; k++) {
                    previ = i - hmm2Parent[k];
                    prevj = j - hmm2Child[k];
                    if (previ >= 0 && prevj >= 0) {
                        if (hmm2Parent[k] != 0) {                // parent present: substitution (* *) or deletion (* -)
                            b = (this == parent.left) ? p.right : p.left;
                            Utils.calcFelsen(felsen, hmm2Child[k] != 0 ? c.seq : null, charTransMatrix, b != null ? b.seq : null, brother.charTransMatrix);
                            if (parent == owner.root || p.orphan || !Utils.USE_UPPER) {
                                emissionProb = Utils.calcEmProb(felsen, equDist);
                            }
                            else {
                            	emissionProb = Utils.calcEmProb(felsen, p.upp);
                            	// If a double deletion (childless parent) then the above
                            	// should just yield the sum of p.upp
                            }
                        } else {                    // insertion (- *)
                            emissionProb = Utils.calcEmProb(c.seq, equDist);
                        }
                        if (previ == 0 && prevj == 0) tr = hmm2PropTransMatrix[START][k];
                        else {
                            for (tr = Utils.log0, prevk = START + 1; prevk < END; prevk++) {
                                tr = Utils.logAdd(tr, probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k]);
                            }
                        }
                        probMatrix[i][j][k] = 1.0*(Math.log(emissionProb) + tr);
                    } else probMatrix[i][j][k] = Utils.log0;
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
    	return hmm2BackProp(1.0);
    }
    double hmm2BackProp(double heat) {
        int emitPatt2State[] = owner.hmm2.getEmitPatt2State();
        final int START = owner.hmm2.getStart();
        final int END = owner.hmm2.getEnd();

        AlignColumn c, p;
        double retVal = 0.0;
        int k, previ, prevj, prevk;
        int pattCode;                                    // alignment column pattern's binary code

        double probMatrix[][][] = hmm2ProbMatrix(heat);

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
    	return hmm2Align(1.0);
    }
    double hmm2Align(double heat) {
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
                prJump[prevk] = heat*(probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k]);
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

        if (Utils.DEBUG) {
        	printToScreenAlignment(0,0);
        }
        return retVal.value;
    }

    /**
     * Samples a new alignment between `this' &amp; `this.parent', taking window sizes into account.
     * Updates all Vertex-stored suppl. data, saving old data for restoration into oldVertex.
     * @return log of 1/proposal
     */
    public double hmm2AlignWithSave() {
    	return hmm2AlignWithSave(1.0);
    }
    public double hmm2AlignWithSave(double heat) {
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
                prJump[prevk] = heat*(probMatrix[previ][prevj][prevk] + hmm2PropTransMatrix[prevk][k]);
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

        if (Utils.USE_UPPER) { 
        	// Then we may have modified the Felsenstein vectors
        	// for `this', so need to restore them.
        	calcFelsen();
        }
        calcOrphan();
        parent.calcFelsen();
        parent.calcOrphan();
        parent.calcIndelLogLike();

//        if (Utils.DEBUG) {
//        	printToScreenAlignment(0,0);
//        }
        
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
        //printToScreenAlignment(0,0,true);
    	StringBuffer sb = new StringBuffer();
    	
    	 owner.root.calcFelsenRecursively();
         owner.root.calcOrphanRecursively();
         owner.root.calcIndelLogLikeRecursively();
         if (Utils.USE_UPPER) {
         	//owner.root.calcFelsenRecursively();
         	owner.root.calcUpperRecursively();
         }   
         //owner.root.recomputeCheckLogLike();
        
        // if (Utils.USE_UPPER) owner.checkUppFelsProducts();
         
    	System.out.println("root.indelLogLike\t "+owner.root.indelLogLike);
        System.out.println("root.orphanLogLike\t "+owner.root.orphanLogLike);

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
        if (Utils.USE_FULL_WINDOWS) winLength = length; p.value = 1.0; 
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
        
        if (Utils.USE_UPPER) {
        	//owner.root.calcFelsenRecursively();
        	owner.root.calcUpperRecursively();
        }        

        double bpp = 0;
        if (!(Utils.USE_FULL_WINDOWS)) {
        	 bpp = -Math.log(p.value / (length - winLength + 1));
        }

        System.out.println("bpp after window select\t "+bpp);

        // select window down (recursively)
        if (left != null && right != null) {
        	left.selectWindow();
        	right.selectWindow();
        }
        selectWindowUp();
	
        owner.root.calcFelsenRecursively();
        owner.root.calcOrphanRecursively();
        owner.root.calcIndelLogLikeRecursively();
        if (Utils.USE_UPPER) {
        	//owner.root.calcFelsenRecursively();
        	owner.root.calcUpperRecursively();
        }   
        
        //printToScreenAlignment(b,b+winLength);
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
//            calcOrphan();
//            parent.calcAllUp();
        } else {
            calcOrphan();
        }
        
        owner.root.calcFelsenRecursively();
        owner.root.calcOrphanRecursively();
        owner.root.calcIndelLogLikeRecursively();
        if (Utils.USE_UPPER) {
        	//owner.root.calcFelsenRecursively();
        	owner.root.calcUpperRecursively();
        }   
        //owner.root.recomputeCheckLogLike();
        //if (Utils.USE_UPPER) owner.checkUppFelsProducts();

        printToScreenAlignment(0,0,true);

        bpp += bppProp;

        System.out.println("root.indelLogLike\t "+owner.root.indelLogLike);
        System.out.println("root.orphanLogLike\t "+owner.root.orphanLogLike);

        if(Utils.DEBUG) {
        	// check proposal - backproposal consistency
        	double bppBack = doRecBackprop();
        	if(parent != null)
        		bppBack += hmm2BackProp();
        	if(Math.abs(bppProp+bppBack) > 1e-5) {
        		System.out.println("Proposal - backproposal inconsistent! Prop: "+bppProp+" Back: "+bppBack);
        	}
        }

        System.out.println("Final window length\t "+winLength);
        double windowProb = 0;
        if (!(Utils.USE_FULL_WINDOWS)) {
        	windowProb = Math.log(Utils.linearizerWeightProb(length, winLength, Utils.WINDOW_MULTIPLIER*Math.sqrt(length)) 
           		/ (length - winLength + 1));
        }
        
        System.out.println("final window selection\t "+windowProb);
        bpp += windowProb;

    	//System.out.println("Total bpp\t "+bpp);
    	//printToScreenAlignment(b,b+winLength);
    	
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
    	owner.root.calcUpperRecursively();
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
     
        //indelLogLike = old.indelLogLike;
//        calcAllUp();
//        owner.root.calcFelsenRecursively();
//        owner.root.calcOrphanRecursively();
//        owner.root.calcIndelLogLikeRecursively();
//        if (Utils.USE_UPPER) {
//        	//owner.root.calcFelsenRecursively();
//        	owner.root.calcUpperRecursively();
//        }   
//        owner.root.recomputeCheckLogLike();
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
    
    /*
    public double fixedColumnNephewUncleSwap() {
    	
    	double logProposalRatio = 0.0;  
    	
        final boolean[] T = {true,false};
    	final int START = owner.hmm2.getStart();
    	final int END = owner.hmm2.getEnd();
    	final int emitPatt2State[] = owner.hmm2.getEmitPatt2State();
       	final int hmm2Parent[] = owner.hmm2.getStateEmit()[0];
        final int hmm2Child[] = owner.hmm2.getStateEmit()[1];
    	
    	final int[][] perms = {{0,1,2},{0,2,1},{2,1,0}};
    	
    	Vertex brother = brother(), uncle = parent.brother(), grandpa = parent.parent, greatgrandpa = grandpa.parent;
    	
    	boolean isLeft = (this == parent.left);
    	boolean uncleIsLeft = (uncle == grandpa.left);
    	boolean gIsRoot = (grandpa.parent == null);
    	//int nTopo = gIsRoot ? 3 : 12;
    	int nTopo = 3;
    	    	
    	calcFelsen();
    	brother.calcFelsen();
    	uncle.calcFelsen();    	
    	if (!gIsRoot) grandpa.calcUpper();
    	
    	ArrayList<double[][]> logMarg = null; // Dynamic programming table
    	double[] logMargTopo = new double[nTopo]; // log marginal probability for each topology
    	
    	ArrayList<double[][]> trans = new ArrayList<double[][]>(3);
    	trans.add(charTransMatrix);
    	trans.add(brother.charTransMatrix);
    	trans.add(uncle.charTransMatrix);
    	
    	boolean[] neighb = new boolean[4]; 
    	// Vector indicating presence/absence of neighbouring nodes
    	// i.e. t, b, u and gg (if the latter exists)
    	
    	// Indexing: t=0, b=1, u=2, gg=3
    	
    	ArrayList<double[]> partialLike = null;
    	
		AlignColumn p = parent.first;
    	AlignColumn g = grandpa.first;
    	
    	double[] fels = new double[owner.substitutionModel.e.length];
    	double[] pUpp = new double[owner.substitutionModel.e.length];
    	
    	AlignColumn t=null, b=null, u=null, gg=null;
    	int columnIndex = 0;
    	int prevk;
        while (p != last || g != grandpa.last) { 
        	// THIS IS NOT CORRECT -- misses out columns that begin with double
        	// gaps in the middle.
        	// -->> would be much easier to first create an array representation of the MSA

        	logMarg.add(new double[nTopo][3]); // Add a new column
        	
        	partialLike = new ArrayList<double[]>(3);
        	
        	// px: 'parent character exists' gx:  'grandpa character exists'
        	// tx: 'this character exists'   bx:  'brother character exists'
        	// ux: 'uncle character exists'  ggx: 'greatgrandpa character exists'
        	boolean px = false, gx = false;
        	boolean tx = false, bx = false; 
        	boolean ux = false, ggx = false;
        	
        	
        	// NEED some more complicated conditions here to ensure that 
        	// px and gx don't become true until we've considered all the intermediate
        	// characters at the tips.
        	if (p.parent != g) { // deletion (* -)
                gx = true;
            } 
        	else if (p.orphan) { // insertion (- *)
                px = true;
            } 
        	else { // substitution (* *)                
                px = true; gx = true;
            }
        	        	
        	if (px) {
        		t = isLeft ? p.left : p.right;
        		tx = (t!=null); 
        		if (tx) {
        			Utils.calcFelsen(fels, t.seq, charTransMatrix, null, null);
        			partialLike.set(0,t.seq);
        		}
        		b = isLeft ? p.right : p.left;
        		bx = (b!=null);
        		if (bx) {
        			Utils.calcFelsen(fels, b.seq, brother.charTransMatrix, null, null);
        			partialLike.set(1,b.seq);
        		}
        	}
        	if (gx) {
    			u = uncleIsLeft ? g.left : g.right;
    			ux = (u!=null); 
    			if (ux) {    				
    				Utils.calcFelsen(fels, u.seq, uncle.charTransMatrix, null, null);
    				partialLike.set(2,fels);
    			}
    			if (!gIsRoot) {
    				ggx = !g.orphan; 
    				if (ggx) {
    					gg = g.parent;    					
    				}
    			}               
    		}        	
        	
        	for (int topo=0; topo<nTopo; topo++) {
	        	int[] perm = perms[topo]; 	        		
	        	neighb[perm[0]] = tx;
	        	neighb[perm[1]] = bx;
	        	neighb[perm[2]] = ux;
	        	if (gIsRoot) {
	        		neighb[3] = false;
	        	}
	        	else {
	        		neighb[3] = ggx;      		
	        		//neighb[perm[3]] = ggx;
	        	}
	        	
	        	// Compute p.upp, to be used in the cases where gx = true
	        	pUpp = new double[owner.substitutionModel.e.length];
	     		for (int i=0; i<pUpp.length; i++) {	     			
	     			for (int j=0; j<pUpp.length; j++) {
	     				pUpp[j] +=  partialLike.get(perm[2])[i] * g.upp[i] * parent.charTransMatrix[i][j];	     				
		        	}
	     		}	        
	        	
    			for (boolean pxx : T) { for (boolean gxx : T) {
    				int pattern = (pxx?1:0) + 2*(gxx?1:0); // (- -)=0, (* -)=1, (- *)=2, (* *)=3    				
    				int state = emitPatt2State[pattern]; 
    				if ((!pxx|px) && (!gxx|gx) && Utils.isValidHistory(pxx,gxx,neighb,gIsRoot)) {    		
    					double logSum = 0;
        				double logEmissionProb = 0.0;
    		        	
        				if (columnIndex==0) ;// need more complicated test for START condition  /////
    	                if (previ == 0 && prevj == 0) logSum = hmm2TransMatrix[START][state];   /////
    	                
    	                else { 
    	                	
    	                	// Sum over possible predecessor indel pair states
    			        	for (logSum = Utils.log0, prevk = START; prevk < END; prevk++) {
    		                    
    			        		double logProb = logMarg.get(columnIndex-1)[topo][prevk];
    			        		
    			        		if (pxx|gxx) {
    			        			// Transition probability of p--g edge  
    			        			logMarg.get(columnIndex-1)[topo][START] hmm2TransMatrix[][state];
    			        			
    			        		    			        		    		                    			        		
    			        		int tState = (neighb[perm[0]]?1:0) + 2*(pxx?1:0);
    			        		int tStatePrev = (neighb[perm[0]]?1:0) + 2*hmm2Child[prevk];
    			        		int bState = (neighb[perm[1]]?1:0) + 2*(pxx?1:0);
    			        		int bStatePrev = (neighb[perm[1]]?1:0) + 2*hmm2Child[prevk];
    			        		
    			        		int uState = (neighb[perm[0]]?1:0) + 2*(gxx?1:0);
    			        		int uStatePrev = (neighb[perm[0]]?1:0) + 2*hmm2Parent[prevk];    			        		
    			        		
    			        		// Add on log probabilities of transitions of the p--t, p--b and g--u branches
    			        		logProb += hmm2TransMatrix[tStatePrev][tState]
    			        		              + hmm2TransMatrix[bStatePrev][bState]
    			        		                 + hmm2TransMatrix[uStatePrev][uState];
    			        		
    			        		if (!gIsRoot) {
    			        			int ggState = (neighb[perm[1]]?1:0) + 2*(gxx?1:0);
    			        			int ggStatePrev = (neighb[perm[1]]?1:0) + 2*hmm2Parent[prevk];
    			        			logProb += hmm2TransMatrix[ggStatePrev][ggState];
    			        		}
    			        		
    			        		logSum = Utils.logAdd(logSum,logProb);        						
    			        	}
    	                }
    		        	
    					if (pxx) { // then either a substitution or an insertion	
    						double[] em = partialLike.get(perm[0]).clone(), em2 = partialLike.get(perm[1]);
    						for (int i=0; i<em.length; i++) em[i] *= em2[i];
    						if (!gxx) { // then an insertion
    							logEmissionProb += Math.log(Utils.calcEmProb(em,owner.substitutionModel.e));
    						}
    						else { // substitution
    							logEmissionProb += Math.log(Utils.calcEmProb(em,pUpp));
    						}    						
    					}
    					else { // Deletion from g to p
    						double tmp = 0.0;
    						for (int i=0; i<g.upp.length; i++) tmp += g.upp[i]; 
    						logEmissionProb = Utils.logAdd(logEmissionProb,Math.log(tmp));
    					}
    		        	logMarg.get(columnIndex)[topo][state] = logEmissionProb + logSum;
    				}
    				else {
    					logMarg.get(columnIndex)[topo][state] = Utils.log0;
    				}
    			} } // End of double loop over possible internal indel states     			    			    		
    		} // End of loop over topologies 
        	
        	if (px) p = p.next;
			if (gx) g = g.next;
			++columnIndex;
        } // End of loop over `columns' 
        
        
        // Compute marginals for each topology, summed over internal indel states 
    	for (int topo=0; topo<nTopo; topo++) {
        	for (int k = START + 1; k < END; k++) {    	
        		logMarg.get(columnIndex-1)[topo][k] += hmm2TransMatrix[k][END];
        		logMargTopo[topo] = Utils.logAdd(logMargTopo[topo],logMarg.get(columnIndex-1)[topo][k]);
        	}
    	}

    	// Select one topology according to its marginal 
    	MuDouble logP = new MuDouble(0.0);
    	int selectedTopo = Utils.logWeightedChoose(logMargTopo, logP);
    	
    	if (selectedTopo == 0) ; // No change to topology
    	else if (selectedTopo == 1) { 
    		/* New topology is
    		        /
    		       g
    		      / \
    		     p   b
    		    / \
    		   t   u
    		   
    		 * /

    		// Set brother and uncle to new values
    		parentNewChild(uncle);
    		parent.parentNewChild(brother);
    	}
    	else if (selectedTopo == 2) {
    		/* New topology is
	                /
	               g
	              / \
	             p   t
	            / \
	           u   b
	   
	         * /
    		
    		// Set this and uncle to new values
    		parentNewChild(uncle);
    		parent.parentNewChild(this);
    	}
    	
    	
    	// Stochastic backtrack to generate new alignment(s) 
    	for (int k = END; k != START; k = prevk) {
            for (prevk = START; prevk < END; prevk++)
                prJump[prevk] = logMarg.get(columnIndex-1)[topo][prevk] + hmm2TransMatrix[prevk][k];
            prevk = Utils.logWeightedChoose(prJump, retVal);
            
            
            t = last;
            b = brother.last;
            p = parent.last;
            g = grandpa.last;
            u = uncle.last;
            
            if (hmm2Parent[prevk] != 0) { // then we need a parent column 
            	p = new AlignColumn(p, true);
            	if () p.left=t;
            	if (neighb[perms[selectedTopo][0]]) p.left=t;
            }
            if (hmm2Child[prevk] != 0) {
            	t = t.prev;
            	b = b.prev;
            }

            p = new AlignColumn(p, true);
            p.seq = 
 
    	    	
    	return logProposalRatio;
    } 
    */
    
    public void nephewUncleRestoreFixedColumns() {
    	
    	Vertex brother = brother(), uncle = parent.brother(), grandpa = parent.parent, greatgrandpa = grandpa.parent;
    	boolean gIsRoot = (grandpa==owner.root);
    	boolean isLeft = (this==parent.left);
    	boolean uncleIsLeft = (uncle==grandpa.left);
         
	    first = old.first;
    	last = old.last;		
    	orphanLogLike = old.orphanLogLike;
        indelLogLike = old.indelLogLike;  
        
        parent.first = parent.old.first;
		parent.last = parent.old.last;
		parent.orphanLogLike = parent.old.orphanLogLike;
		parent.indelLogLike = parent.old.indelLogLike;
        
		brother.first = brother.old.first;
		brother.last = brother.old.last;
		brother.orphanLogLike = brother.old.orphanLogLike;
		brother.indelLogLike = brother.old.indelLogLike;
        
		grandpa.first = grandpa.old.first;
		grandpa.last = grandpa.old.last;
		grandpa.orphanLogLike = grandpa.old.orphanLogLike;
		grandpa.indelLogLike = grandpa.old.indelLogLike;
        
		uncle.first = uncle.old.first;
		uncle.last = uncle.old.last;
		uncle.orphanLogLike = uncle.old.orphanLogLike;
		uncle.indelLogLike = uncle.old.indelLogLike;
		
		if (!gIsRoot) {
			greatgrandpa.first = greatgrandpa.old.first;
			greatgrandpa.last = greatgrandpa.old.last;
			greatgrandpa.orphanLogLike = greatgrandpa.old.orphanLogLike;
			greatgrandpa.indelLogLike = greatgrandpa.old.indelLogLike;
		}
		
		uncle.parent = parent; 		
 		//uncle.last.parent = parent.last;
 		if (isLeft)		 { parent.left = uncle; }//parent.last.left = uncle.last; } 
 		else 			 { parent.right = uncle; }//parent.last.right = uncle.last; }
        parent = grandpa;
 		//last.parent = grandpa.last;		 				 		
 		if (uncleIsLeft) { grandpa.left = this; }//grandpa.last.left = last; } 
 		else 		 	 { grandpa.right = this; }//grandpa.last.right = last; }
 		
 		
 		if (Utils.DEBUG) {        	
 			uncle.printPointers(); //uncle.printPointers2();
         	brother.printPointers();  //brother.printPointers2();         	
         	if (!gIsRoot) { uncle.parent.printPointers(); }//parent.printPointers2(); }
         	printPointers();
         	uncle.parent.printPointers();//uncle.parent.printPointers2();
         	
         	         	  
         	//owner.checkPointers();
         }
         
         
        DEBUG=2; // Activate verbose print statements 
     	uncle.calcAllUp(); // uncle is now lower 
     	DEBUG=0; // Deactivate verbose print statements
    	
    }
    
    public double nephewUncleSwapFixedColumns() {
    	double logProposalRatio = 0.0;
    	
    	Vertex brother = brother(), uncle = parent.brother(), grandpa = parent.parent, greatgrandpa = grandpa.parent;
    	boolean gIsRoot = (grandpa==owner.root);
    	boolean isLeft = (this==parent.left);
    	boolean uncleIsLeft = (uncle==grandpa.left);
    	boolean grandpaIsLeft = false; 
    	if (!gIsRoot) grandpaIsLeft = (grandpa==greatgrandpa.left);

        
    	  if (Utils.DEBUG) {
          	printPointers(); //printPointers2();        	
          	brother.printPointers();  //brother.printPointers2();
          	if (!gIsRoot) { parent.printPointers(); }//parent.printPointers2(); }
          	uncle.printPointers(); //uncle.printPointers2();
          	parent.printPointers();//parent.printPointers2();  
          	//owner.checkPointers();
          }
    	
    	// Extract full alignment for all nodes.
    	// NB this also has the effect of saving the old state, which is useful, albeit slow.
    	String[] ali = owner.getState().getFullAlign();
    	
    	double P = 10.2; // Probability of an additional indel being incorporated
    	// As it stands, P is not really a probability. Also, should rewrite the 
    	// references to this below in logProposalRatio so that we can have
    	// more general forms for the weights.
    	double P2 = Math.pow(P,2);
    	// Should really derive this from lambda, mu and R somehow
		double[] weights2 = {1,P};
		double[] weights3 = {1,P,P2};
        		
		// Begin saving of current Vertex objects
		old.last = last.clone();
		old.orphanLogLike = orphanLogLike;
        old.indelLogLike = indelLogLike;
        old.parent = parent.old; // for printing
		
        brother.old.last = brother.last.clone();
        brother.old.orphanLogLike = brother.orphanLogLike;
        brother.old.indelLogLike = brother.indelLogLike;
        brother.old.parent = parent.old; // for printing
        
        parent.old.last = parent.last.clone();
        old.last.parent = parent.old.last;
        brother.old.last.parent = parent.old.last;
        if (isLeft) { parent.old.last.left = old.last; parent.old.last.right = brother.old.last; }
        else 	    { parent.old.last.right = old.last; parent.old.last.left = brother.old.last; }
		parent.old.orphanLogLike = parent.orphanLogLike;
		parent.old.indelLogLike = parent.indelLogLike;
		parent.old.parent = grandpa.old; //	for printing
		
		uncle.old.last = uncle.last.clone();
        uncle.old.orphanLogLike = uncle.orphanLogLike;
        uncle.old.indelLogLike = uncle.indelLogLike;
        uncle.old.parent = grandpa.old; // for printing
        
        grandpa.old.last = grandpa.last.clone();
        parent.old.last.parent = grandpa.old.last;
        uncle.old.last.parent = grandpa.old.last;
        if (uncleIsLeft) { grandpa.old.last.left = uncle.old.last; grandpa.old.last.right = parent.old.last; }
        else 	         { grandpa.old.last.right = uncle.old.last; grandpa.old.last.left = parent.old.last; }
        grandpa.old.orphanLogLike = grandpa.orphanLogLike;
        grandpa.old.indelLogLike = grandpa.indelLogLike;        
		        
        
        if (!gIsRoot) {
            greatgrandpa.old.last = greatgrandpa.last.clone();
            grandpa.old.last.parent = greatgrandpa.old.last;
            if (grandpaIsLeft) greatgrandpa.old.last.left = grandpa.old.last; 
            else 	           greatgrandpa.old.last.right = grandpa.old.last;
        	greatgrandpa.old.orphanLogLike = greatgrandpa.orphanLogLike;
        	greatgrandpa.old.indelLogLike = greatgrandpa.indelLogLike;
        	grandpa.old.parent = greatgrandpa.old; // for printing
        }
   		
        // Loop over columns of the initial alignment, starting from the last (non-virtual) column
        AlignColumn t=last.prev, p=parent.last.prev, b=brother.last.prev, g=grandpa.last.prev, u=uncle.last.prev;
        AlignColumn gg = null; if (!gIsRoot) gg=greatgrandpa.last;
        
        // Old columns, for restoration
        AlignColumn to=old.last, po=parent.old.last, bo=brother.old.last, go=grandpa.old.last, uo=uncle.old.last;   
        AlignColumn ggo = null; if (!gIsRoot) ggo=greatgrandpa.old.last;
 		
        for (int col=ali[0].length()-1; col>=0; col--) {

        	System.out.println();
        	System.out.print(ali[brother.index].charAt(col));        	
        	System.out.print(ali[index].charAt(col));
        	System.out.print(ali[uncle.index].charAt(col));
        	System.out.print(ali[parent.index].charAt(col));
        	System.out.print(ali[grandpa.index].charAt(col));
        	System.out.println();
        	boolean tx = (ali[index].charAt(col)!='-'); 		
        	boolean px = (ali[parent.index].charAt(col)!='-'); 	
        	boolean bx = (ali[brother.index].charAt(col)!='-'); 
        	boolean gx = (ali[grandpa.index].charAt(col)!='-');
        	System.out.println(gx+" "+ali[grandpa.index].charAt(col));
        	boolean ux = (ali[uncle.index].charAt(col)!='-');
          	boolean ggx = false;
        	if (!gIsRoot) ggx = (ali[greatgrandpa.index].charAt(col)!='-');
        	        	
        	int nTips = (tx?1:0) + (bx?1:0) + (ux?1:0) + (ggx?1:0);
        	System.out.println("col = "+col+", nTips = "+nTips);
        	
        	System.out.println("g = "+g.toString());
//        	if (Utils.DEBUG) {
//        		System.out.println(ali[index].charAt(col));
//        		System.out.println(ali[parent.index].charAt(col));
//        		System.out.println(ali[brother.index].charAt(col));
//        		System.out.println(ali[grandpa.index].charAt(col));
//        		System.out.println(ali[uncle.index].charAt(col));
//        	}        	      
        	
        	// Save current columns      	        	        	       	
        	if (!gIsRoot && ggx) {
        		ggo.prev = gg.clone(); 
        		ggo.prev.next = ggo; gg = ggo.prev;        		    			    			
        	}
        	if (gx) {        		
        		go.prev = g.clone(); 
        		go.prev.next = go; go = go.prev;
        		go.parent = gIsRoot ? null : ggo;
        		if (!gIsRoot && ggx) {        			
    				if (grandpaIsLeft)  ggo.left = go;
    				else 		 	    ggo.left = go;
    			}        		
        	}
        	if (ux) { 
        		uo.prev = u.clone(); 
        		uo.prev.next = uo; uo = uo.prev;
        		uo.parent = go;
        		if (gx) {        				
        			if (uncleIsLeft)  go.left = uo;
    				else 		 	  go.right = uo;
        		}
        	}
        	if (px) {
        		po.prev = p.clone(); 
        		po.prev.next = po; po = po.prev;
        		po.parent = go;
        		if (gx) {        			
        			if (uncleIsLeft)  go.right = po;
    				else 		 	  go.left = po;
        		}        		
        	}
        	if (bx) { 
        		bo.prev = b.clone(); 
        		bo.prev.next = bo; bo = bo.prev;
				bo.parent = po;
        		if (px) {    			
    				if (isLeft)  po.right = bo;
    				else 		 po.left = bo;
    			}           		
        	}
        	if (tx) { 
        		to.prev = t.clone(); 
        		to.prev.next = to; to = to.prev;
        		to.parent = po;
        		if (px) {        			
    				if (isLeft)  po.left = to;
    				else 		 po.right = to;
    			}        		
        	}
        	
        	if (Utils.DEBUG) {
	        	if (px&gx) if (p==g) System.out.println("BEFORE: p == g? Should be false. Actually is "+(p==g));
		    	if (tx&gx) {
		    		if (px) if (t.parent!=p) System.out.println("BEFORE: t's parent is p? Should be true. Actually is "+(t.parent==p));
		    	   	if (t.parent==g) System.out.println("BEFORE: t's parent is g? Should be false. Actually is "+(t.parent==g));
		    	}
		    	if (ux&px) {
		    		if (gx) if (u.parent!=g) System.out.println("BEFORE: u's parent is g? Should be true. Actually is "+(u.parent==g));
		    		if (u.parent==p) System.out.println("BEFORE: u's parent is p? Should be false. Actually is "+(u.parent==p));
		    	}	   
	        }
	    	
        	// Consider the possible cases
        	if (nTips == 3 || nTips == 4 || ( nTips == 2 && ((bx&ggx)|(tx&ux)) ) ) {        	
        		// Then we have a tip present on both sides of the central edge,
        		// both before and after the move, so we need not do anything
        		// except for changing the parentage (and recomputing the 
        		// emission probabilities).        	
        	}        	
        	else if ((tx&ggx)|(bx&ux)) {
    			// Case 1:
    		    //    grandpa exists before and after
    			//    parent exists before, but possibly not after
        		//    u does not exist
        		// Case 2:
        		//    parent exists before and after
    			//    grandpa exists before, but possibly not after
        		//    t does not exist
        		
        		//    log(bwd) = 0
    			      		
        		logProposalRatio += Math.log(1+P); // denominator for -log(fwd)
        		
    			if (Utils.weightedChoose(weights2)==0) {    				
    				if (tx) { delete(p); p = p.prev; px = false; }
    				else    { delete(g); g = g.prev; gx = false; }
    				// log(fwd) = log(1)
    			}
    			else logProposalRatio -= Math.log(P);    			   			
        	}        	
        	else if ((bx&tx)|(ux&ggx)) {
        		// Case 1:
    		    //    parent exists before and after
    			//    grandpa exists after, and possibly before
        		//    u and gg do not exist
        		// Case 2:
    		    //    grandpa exists before and after
    			//    parent exists after, and possibly before
        		//    b and t do not exist

        		// -log(fwd) = 0
        		        		
        		logProposalRatio -= Math.log(1+P); // denominator of log(bwd)
        		
        		if (bx) {
        			if (gx) logProposalRatio += Math.log(P); // grandpa existed before
        			else { // Insert new g column        			
            			g = new AlignColumn(g.next); // New orphan column          			
            			gx = true;
            			// t is already not marked as orphan because p exists
            		}
        		}
        		else {
        			if (px) { // If parent existed before
            			logProposalRatio += Math.log(P);
            		}
            		else { // Insert new p column
            			p = new AlignColumn(p.next); 
            			p.parent = g;
            			p.orphan = false;
            			px = true;
            		}	
        		}
        		
        	}        
        	else if (nTips == 2) {
        		throw new RuntimeException("We missed a case where nTips==2.");
        	}
       	
        	// Now on to the cases where there was only one tip character to 
        	// begin with:
        	
        	else if (bx|ggx) {
        		System.out.println("bx|ggx");
        		// Then the only character in this column was at the brother
        		// or the greatgrandpa, in which case there are no nephew or
        		// uncle characters to swap, so there is nothing to do.  
        		// We could sample new internal states, but it's not necessary,
        		// so we will leave any such columns as they are.
        	}
        	else if (tx|ux) {
        		System.out.println("tx|ux");
        		// Then:        		
    		    //    parent exists possibly before and possibly after
    			//    grandpa exists possibly before and possibly after
        		
        		// Back proposal probabilities
        		if (px&gx)  logProposalRatio -= 2*Math.log(P);
    			if (px&!gx) logProposalRatio -= Math.log(P);
    			if (!px); //logProposalRatio -= Math.log(1);
    				
        		// First sample new internal states        		
        		int choice = Utils.weightedChoose(weights3);
        		if (choice == 0) { // Then we're choosing no internals
    			    // -log(fwd) = log(1) = 0
        			
        			if (px) { delete(p); p = p.prev; px = false; }  
        			if (gx) { delete(g); g = g.prev; gx = false; } 
        			// I don't think the order of the above steps matters
        			
        			if (tx) t.orphan = true;
    				if (ux) u.orphan = true;
        		}
        		else if (choice == 1) { // Then we're choosing to impute only one character    				
    				// -log(fwd) = log(P)
    				logProposalRatio += Math.log(P);
    				
    				if (tx) {
    					// if we're considering t, then its new parent must be g
    					if (px) { delete(p); p = p.prev; px = false; } 
	    				if (!gx) {
	    					g = new AlignColumn(g.next); // New orphan column   
	    					gx = true;
	    				}    				
	    				t.orphan = false;
    				}
    				else if (ux) {
    					// if we're considering u, then its new parent must be p
    					if (gx) { delete(g); g = g.prev; gx = false; }
	    				if (!px) {
	    					p = new AlignColumn(p.next); // New orphan column   
	    					px = true;
	    				}    				
	    				u.orphan = false;
    				}
    				
    			}
        		else if (choice == 2) { //Then we're imputing characters at p and g
    				// -log(fwd) = 2 log(P) 
    				logProposalRatio += 2*Math.log(P);
    				
    				if (!gx) {
    					System.out.println("g.prev = "+g.prev.toString());
    					g = new AlignColumn(g.next); // New orphan column
    					System.out.println("g.prev = "+g.prev.toString());
    					gx = true;
    				}
    				if (!px) {
    					p = new AlignColumn(p.next);            			
            			p.orphan = false;
            			px = true;
    				}    				
    				if (tx) t.orphan = false;
    				if (ux) u.orphan = false;
    			}
        	}        	
        	else {
        		System.out.println("None of the above");
        		throw new RuntimeException("We reached a case that wasn't handled properly.");
        	}
        	
	    	if (col==0) {
	    		first = t; old.first = to;         		 
	       		parent.first = p; parent.old.first = po; 
	      		brother.first = b; brother.old.first = bo; 
	      		uncle.first = u; uncle.old.first = uo; 
	      		grandpa.first = g; grandpa.old.first = go;
	      		System.out.println("first = "+first.toString()+", t = "+t.toString());
	      		System.out.println("parent.first = "+parent.first.toString()+", p = "+p.toString());
	      		System.out.println("brother.first = "+brother.first.toString()+", b = "+b.toString());
	      		System.out.println("uncle.first = "+uncle.first.toString()+", u = "+u.toString());
	      		System.out.println("grandpa.first = "+grandpa.first.toString()+", g = "+g.toString());
          	}
        	  
	    	
        	// Reassign any parent pointers. 
	    	// If the character is to be an orphan, we specify the next column
	    	// in the parent sequence as the adopted parent.
	   if (tx) System.out.println("t = "+t.toString()+", t.parent = "+t.parent.toString()+", t.parent.left = "+t.parent.left.toString());
        	if (tx) t.updateParent(t.orphan?g.next:g,uncleIsLeft);
       if (tx) System.out.println("t = "+t.toString()+", t.parent = "+t.parent.toString()+", t.parent.left = "+t.parent.left.toString());
        	if (px) p.updateParent(p.orphan?g.next:g,!uncleIsLeft);
       if (tx) System.out.println("t = "+t.toString()+", t.parent = "+t.parent.toString()+", t.parent.left = "+t.parent.left.toString());
        	if (ux) u.updateParent(u.orphan?p.next:p,isLeft);
       if (tx) System.out.println("t = "+t.toString()+", t.parent = "+t.parent.toString()+", t.parent.left = "+t.parent.left.toString());
        	if (bx) b.updateParent(b.orphan?p.next:p,!isLeft);
       if (tx) System.out.println("t = "+t.toString()+", t.parent = "+t.parent.toString()+", t.parent.left = "+t.parent.left.toString());
        	if (gx & !gIsRoot) g.updateParent(g.orphan?gg.next:gg,grandpaIsLeft);        
        	
        	if (Utils.DEBUG) {
	        	if (px&gx) if (p==g) System.out.println("AFTER: p == g? Should be false. Actually is "+(p==g));
		    	if (tx&gx) {
		    		if (px) if (t.parent==p) System.out.println("AFTER: t's parent is p? Should be false. Actually is "+(t.parent==p));
		    	   	if (t.parent!=g) System.out.println("AFTER: t's parent is g? Should be true. Actually is "+(t.parent==g));
		    	}
		    	if (ux&px) {
		    		if (gx) if (u.parent==g) System.out.println("AFTER: u's parent is g? Should be false. Actually is "+(u.parent==g));
		    		if (u.parent!=p) {
		    			System.out.println("AFTER: u's parent is p? Should be true. Actually is "+(u.parent==p));
		    			System.out.println(u.toString()+" "+u.parent.toString()+" "+p.toString());
		    		}
		    	}	   
	        }
        	System.out.println("*"+gx+" "+ali[grandpa.index].charAt(col));

        	// Decrement the column pointers
        	if (tx) t = t.prev; 
        	if (px) p = p.prev;
        	if (bx) b = b.prev;
        	if (gx) g = g.prev;        	
        	if (ux) u = u.prev;
        	if (ggx) gg = gg.prev;       	
        	
        }   
        
        // Swap the pointers between vertices, and the first and last columns
        // NB the uncle has to be done first, otherwise parent refers to grandpa
        // NB this also has to be done after all other references to parent in reference to the old parent
        // (or otherwise we could instead refer to a pointer to the parent under a different name,
        // rather than using the member variable `parent')
        uncle.parent = parent; 		
 		uncle.last.parent = parent.last;
 		if (isLeft)		 { parent.left = uncle; parent.last.left = uncle.last; } 
 		else 			 { parent.right = uncle; parent.last.right = uncle.last; }
        parent = grandpa;
 		last.parent = grandpa.last;		 				 		
 		if (uncleIsLeft) { grandpa.left = this; grandpa.last.left = last; } 
 		else 		 	 { grandpa.right = this; grandpa.last.right = last; }
 		
        if (Utils.DEBUG) {
        	old.printPointers();        	
        	brother.old.printPointers(); 
        	if (!gIsRoot) parent.old.printPointers(); 
        	uncle.old.printPointers(); 
        	uncle.parent.old.printPointers(); 
        	System.out.println("------------------");
        	printPointers(); //printPointers2()
        	brother.printPointers();  //brother.printPointers2();
        	if (!gIsRoot) { parent.printPointers(); }//parent.printPointers2(); }
        	uncle.printPointers(); //uncle.printPointers2();
        	uncle.parent.printPointers();//uncle.parent.printPointers2();
        	
        	owner.checkPointers();        	
        }
        
        
        DEBUG=2; // Activate verbose print statements 
    	uncle.calcAllUp(); // uncle is now lower
    	DEBUG=0; // Deactivate verbose print statements
    	
//    	if (gIsRoot) grandpa.calcUpperRecursively();
//    	else         greatgrandpa.calcUpperRecursively();
    	// The recalculation of upp should only be done inside
    	// a function that needs those vectors, since they are
    	// not needed for the computation of the likelihood.    
    		
    	return logProposalRatio;
    }
    /**
     * Swaps {@code this} with its uncle, but proposes new alignments only between this node and its parent, and
     * the uncle node and its parent, every other alignment is kept fixed. Slow because full sequence
     * alignment is done.
     * Assumes {@code this} has a non-null grandparent.
     * @return log-quotient of backproposal and proposal
     */
    public double swapWithUncleAlignToParent() {
       	boolean old_USE_UPPER = Utils.USE_UPPER;
       	int INDEX = -1;
       	if (index == INDEX) {
       		Utils.USE_UPPER = false;
       	}
       	System.out.println("************"+old_USE_UPPER);
    	//Utils.USE_UPPER = true;
       	//Utils.USE_UPPER = false;
       	if (Utils.USE_UPPER) owner.root.calcUpperRecursively();
    	Vertex uncle = parent.brother(), grandpa = parent.parent;

        System.out.println("this = "+index+", uncle = "+uncle.index);
        System.out.println(owner.printedTree());
        
        calcAllUp(); ///
       	owner.root.calcUpperRecursively(); ///
       	

        // make node and window selection
        lastSelected();
        uncle.lastSelected();
        // Joe: why do we need the above, when this function only uses hmm2?
        
        fullWin();
        parent.fullWin();
        uncle.fullWin();
        grandpa.fullWin();
       	
        printToScreenAlignment(0,0,true);
        
        double ret = 0.0;

        // Choose whether to realign uncle or nephew first.
        // We choose according to how large the subtrees are
        // for the nephew and the uncle, because we want to include
        // as much information as possible about the rest of the tree
        // for the second realignment step.
        // If Utils.USE_UPPER is false, then we always choose to realign
        // the lower subtree first, because the alternative would
        // involve no information being used.
		double[] weights = new double[2];
		boolean realignLowerFirst = true;
		if (Utils.USE_UPPER) {
	        parent.calcUpperWithCheck();
			owner.countLeaves();
	        weights[0] = Math.pow(leafCount,Utils.LEAF_COUNT_POW);
	        weights[1] = Math.pow(uncle.leafCount,Utils.LEAF_COUNT_POW);
	       // realignLowerFirst = (Utils.weightedChoose(weights) == 1);
	        
	        // CURRENTLY hmm2AlignWithSave operates in such a way 
	        // that it doesn't work if the upper node is realigned first
	        realignLowerFirst = true;
	        
		    if (realignLowerFirst) {
		    	System.out.println("Aligning lower first.");
		    	parent.calcUpperWithoutBrother();

		    	//ret += Math.log(weights[1]/weights[0]);
		    }
		    else {
		    	System.out.println("Aligning upper first.");
		    	parent.calcFelsenWithout(this);
		    	//ret += Math.log(weights[0]/weights[1]);
		    }
		}
        System.out.println("log move order probability ratio = "+ret);
        
        owner.root.calcIndelLogLikeRecursively();

        // compute alignment backproposal
        double bpp = uncle.hmm2BackProp(PROPOSAL_HEAT);
	    System.out.println("After first back prop = "+bpp);    
        bpp += hmm2BackProp(PROPOSAL_HEAT);
        System.out.println("log back proposal probability = "+bpp);
        ret += bpp;

        // do the swap
        parentNewChild(uncle);			// order is important here
        uncle.parentNewChild(this);
        uncle.parent = parent;
        parent = grandpa;
     

        // align the sequences
        double bppProp = 0.0;
        if (realignLowerFirst) {
        	bppProp += uncle.hmm2AlignWithSave(PROPOSAL_HEAT); // This calls uncle.parent.calcFelsen()
        	//uncle.parent.calcFelsenWithCheck();

        	//bppProp += uncle.hmm2Align(); // This calls parent.calcFelsen()
        	//uncle.parent.calcFelsenWithCheck();
        	System.out.println("After first realignment = "+bppProp);
        	bppProp += hmm2AlignWithSave(PROPOSAL_HEAT);
        	//bppProp += hmm2Align();
        	//parent.calcFelsenWithCheck();
        	//parent.calcFelsen();
        }
        else {
        	bppProp += hmm2AlignWithSave(PROPOSAL_HEAT);
        	if (Utils.USE_UPPER) {
        		uncle.parent.calcUpper();
        	}
        	bppProp += uncle.hmm2AlignWithSave(PROPOSAL_HEAT);
        }
        
        
        System.out.println("log forward proposal probability = "+(-bppProp));
        ret += bppProp;
        
        uncle.calcAllUp();         
        owner.root.calcIndelLogLikeRecursively();

        owner.root.calcUpperRecursively();
        
//        Utils.USE_UPPER = false;
//        double bppProp2 = uncle.hmm2AlignWithSave(); // This calls parent.calcFelsen()
//    	System.out.println("After first realignment = "+bppProp2);
//    	bppProp2 += hmm2AlignWithSave();
//    	System.out.println("log forward proposal probability = "+(-bppProp2));
//    	Utils.USE_UPPER = old_USE_UPPER;
        // After this, all the seq and upp vectors will include
        // all the relevant contributions.

        if(Utils.DEBUG) {
        	System.out.println(owner.printedTree());
        	// check proposal - backproposal consistency
        	if (realignLowerFirst) {
		    	uncle.parent.calcUpperWithoutBrother();
		    }
		    else {
		    	uncle.parent.calcFelsenWithout(uncle);
		    }
        	double bppBack = uncle.hmm2BackProp(PROPOSAL_HEAT);
        	bppBack += hmm2BackProp(PROPOSAL_HEAT);
        	if (realignLowerFirst) {
		    	uncle.parent.calcUpper();
		    }
		    else {
		    	uncle.parent.calcFelsen();
		    }
        	if(Math.abs(bppProp+bppBack) > 1e-8) {
        	  System.out.println("###Proposal - backproposal inconsistent in swapWithUncleAlignToParent! Prop: "+bppProp+" Back: "+bppBack);
        	}
//        	uncle.calcAllUp();
//        	owner.root.calcUpperRecursively();
//        	bppBack = uncle.hmm2BackProp();
//        	bppBack += hmm2BackProp();
//        	
//        	if(Math.abs(bppProp+bppBack) > 1e-8) {
//        	  System.out.println("######Proposal - backproposal inconsistent in swapWithUncleAlignToParent! Prop: "+bppProp+" Back: "+bppBack);
//        	}
        }
        //owner.changingTree = false;
    	Utils.USE_UPPER = old_USE_UPPER;
        printToScreenAlignment(0,0,true);

        if (index == INDEX) {
        	System.out.println("Proposal ratio = "+ret);
        	return Double.POSITIVE_INFINITY;
        }
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
        
        owner.root.calcFelsenRecursively();
        owner.root.calcOrphanRecursively();
        owner.root.calcIndelLogLikeRecursively();
        if (Utils.USE_UPPER) {
        	//owner.root.calcFelsenRecursively();
        	owner.root.calcUpperRecursively();
        }   
       // owner.root.recomputeCheckLogLike();
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
    	boolean old_USE_UPPER = Utils.USE_UPPER; 
    	Utils.USE_UPPER = false;
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

        Utils.USE_UPPER = old_USE_UPPER;
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
            v.calcIndelLogLike();
        }
//        if (Utils.USE_UPPER) {
//        	owner.root.calcUpperRecursively();
//        }
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
    
    public String[] printedAlignment(boolean inWindow) {
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
    public String[] printedAlignment() {
    	return printedAlignment(false);
    }
    String[] printedAlignmentInWindow() {
    	return printedAlignment(true);
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
    
    void printPointers2() {
        AlignColumn c = last;
        AlignColumn p = parent.last;
        while (p != null) {
            System.out.print(p + " ");
            p = p.prev;
        }
        System.out.println();
        while (c != null) {
            System.out.print(c.parent + " ");
            c = c.prev;
        }
        System.out.println("\n");

        c = last.prev;
        p = parent.last.prev;
        while (p.prev != null) {
            System.out.print(p.mostLikely() + "");
            p = p.prev;
        }
        System.out.print(p.mostLikely() + "");
        System.out.println();
        while (c != null) {
            if (c.parent.next != null) {
                System.out.print(c.parent.mostLikely() + "");
            }
            c = c.prev;
        }
        System.out.println();
        c = last.prev;
        while (c != null) {
            if (c.parent.next != null) {
                System.out.print((c.orphan ? "y" : "n"));
            }
            c = c.prev;
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
    	printToScreenAlignment(windowStart, windowEnd,false);
    }
    public void printToScreenAlignment(int windowStart, int windowEnd, boolean printAll) {
		int nGapsInWindow = 0;
    	if (parent != null) {
    		String[] s0 = printedAlignmentInWindow();
    		System.out.println("parent\t"+s0[0]);
    		System.out.println("child\t"+s0[1]+"\n");
//    		for (int k=0; k< windowEnd; k++) {
//    			nGapsInWindow += (s0[1].charAt(k) == '-') ? 1 : 0;
//    		}
    	}
    	if (!printAll) return;
		String[] s = owner.getState().getFullAlign();
		if (windowEnd == 0) {
			for (String sequence : owner.getState().seq) {
				windowEnd = Math.min(sequence.length(),windowEnd);
			}
		}
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
    		String n = owner.vertex[ii].name;
    		if (n == null) {
    			//n = String.format("%-8.4f",owner.vertex[ii].edgeLength);
    			n = String.format("%-8d",owner.vertex[ii].index);
    		}
    		System.out.println(String.format("%-8s", n)+"\t"+s[ii]);
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
                throw new Error("Problem is vertex " + index + ":\np is: " + p + " p.left is " + p.left + " p.left.orphan: " + p.left.orphan +
                        " p.left.parent: " + p.left.parent);
            }
            if (p.right != null && (p.right.orphan || p.right.parent != p)) {
            	System.out.println(p.owner.index+" "+p.prev.toString()+" "+p.prev.mostLikely());
                throw new Error("Problem is vertex " + index + ":\np is: " + p + " p.right is " + p.right + " p.right.orphan: " + p.right.orphan +
                        " p.right.parent: " + p.right.parent);
            }
        }
        for (AlignColumn l = left.first; l != null; l = l.next) {
            if (!l.orphan) {
            	System.out.println("l = "+l.toString()+", p = "+l.parent.toString());
                if (l.parent == null || l.parent.left != l) {
                	System.out.println("An error is about to occur.........");
                	l.owner.printPointers();
                    throw new Error("Problem in vertex "+index+"/"+l.owner.index+ ":\nl is: " + l + " ("+l.mostLikely()+") " + (l.parent == null ? " l does not have a parent" : " l parent is: " + l.parent + " l parent left is: " + l.parent.left));
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

//    public static void main(String[] args){
////    	boolean[] tips = {true,true,false};
////    	boolean[] x = {true,false};
////    	for (boolean p : x) {
////    		for (boolean g : x) {
////    			System.out.print(Utils.isValidTopology(p,g,tips,true)+" ");
////    		}    		
////    		System.out.println();
////    	}
////    	ArrayList<Double> a = new ArrayList<Double>(5);
////    	System.out.println("a.size() = "+a.size());
////    	a.add(0.0); a.add(1.0); a.add(2.0);
////    	ArrayList<ArrayList<Double> > x = new ArrayList<ArrayList<Double> >();
////    	x.add(a);
////    	
////    	ArrayList<Double> b = x.get(0);
////    	b.set(0,238.0);
////    	System.out.println("a = "+a.get(0)+" "+a.get(1)+" "+a.get(2));
////    	System.out.println("b = "+b.get(0)+" "+b.get(1)+" "+b.get(2));
////    	System.out.println("x.get(0) = "+x.get(0).get(0)+" "+x.get(0).get(1)+" "+x.get(0).get(2));
//    	
//        	
//    }
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
