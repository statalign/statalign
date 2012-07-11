package statalign.postprocess.utils;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Michael
 */
public class RNAFoldingTools {
	
    public static void main(String[] args) {

    	// Generate an example "Base pairing probability" matrix
    	String [] structures = {"((((...))))....", 
    							"((((...))))(.).", 
    							"((((...)))).(.)"};
    	double[][] basePairCount = getBasePairCountMatrix(structures);
    	//double[] singleBaseCount = getSingleBaseCount(structures);
    	
    	// Perform the posterior decoding
    	RNAFoldingTools rnaTools = new RNAFoldingTools();
    	int [] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(basePairCount); // fails for this example, because not base-pairing probability matrix
    	//int [] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(basePairCount, singleBaseCount); // this works
    	
    	// Print the list of paired sites for the posterior decoding structure
    	for(int i = 0 ; i < pairedSites.length ; i++)
    	{
    		System.out.println((i+1)+"\t"+pairedSites[i]);
    	}
    	
    	// Print a dot bracket string representation of the posterior decoding structure
    	System.out.println("Dot bracket structure: " + getDotBracketStringFromPairedSites(pairedSites));
    	//System.out.println("Dot bracket structure (sanity check): " + getDotBracketStringFromPairedSites(getPairedSitesFromDotBracketString(getDotBracketStringFromPairedSites(pairedSites))));

    }

    /**
     * General purpose class for holding a pair of integers.
     */
    static class Pair {

        int x;
        int y;

        public Pair(int x, int y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public String toString() {
            return "(" + x + ", " + y + ")";
        }
    }
    public static final double emptyValue = Double.MAX_VALUE;

    /**
     * A single-threaded method for generating the posterior-decoding structure.
     * @param basePairProb a NxN matrix of base-pairing probabilities.
     * @param singleBaseProb an array of length N representing probabilities for unpaired bases.
     * @return an array of paired positions. Where (i, array[i]) represents a nucleotide pairing between nucleotides (i+1, array[i]), if array[i] = 0, then (i+1) is unpaired.
     */
    public static int [] getPosteriorDecodingConsensusStructure(double[][] basePairProb, double[] singleBaseProb) {
        double[][] eMatrix = new double[basePairProb.length][basePairProb[0].length];
        for (int i = 0; i < eMatrix.length; i++) {
            for (int j = 0; j < eMatrix.length; j++) {
                eMatrix[i][j] = RNAFoldingTools.emptyValue;
            }
        }
        
        int[] pairedWith = new int[eMatrix.length];
        recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, 0, eMatrix.length - 1, pairedWith);

        return pairedWith;
    }

    /**
     * Performs posterior-decoding using multi-threading.
     * @param basePairProb  a NxN matrix of base-pairing probabilities.
     * @param singleBaseProb an array of length N representing probabilities for unpaired bases.
     * @return an instance of MultiThreadedPosteriorDecoding, which has public access to various useful variables (e.g. paired nucleotides of the MEA structure).
     */
    public MultiThreadedPosteriorDecoding performPosteriorDecodingMultiThreaded(double[][] basePairProb, double[] singleBaseProb) {
        MultiThreadedPosteriorDecoding m = new MultiThreadedPosteriorDecoding(basePairProb, singleBaseProb);
        return m;
    }
    
    /**
     * Performs posterior-decoding using multi-threading.
     * @param basePairProb a NxN matrix of base-pairing probabilities. Assumes that sum(row) <= 1, in order to calculate the single base probabilities.
     * @return an instance of MultiThreadedPosteriorDecoding, which has public access to various useful variables (e.g. paired nucleotides of the MEA structure).
     */
    public MultiThreadedPosteriorDecoding performPosteriorDecodingMultiThreaded(double[][] basePairProb) {
    	double [] singleBaseProb = new double[basePairProb.length];
    	for(int i = 0 ; i < basePairProb.length ; i++)
    	{
    		singleBaseProb[i] = 1;
    		for(int j = 0 ; j < basePairProb[0].length ; j++)
        	{
    			singleBaseProb[i] -= basePairProb[i][j];
        	}    		
    	}

    	for(int i = 0 ; i < singleBaseProb.length ; i++)
    	{
    		System.out.println((i+1)+"\t"+singleBaseProb[i]);
    	}
    	
    	return new MultiThreadedPosteriorDecoding(basePairProb, singleBaseProb);
    }
    
    /**
     * Returns the posterior-decoding consensus structure.
     * @param basePairProb a NxN matrix of base-pairing probabilities.
     * @return the posterior-decoding consensus structure.
     */
    public int [] getPosteriorDecodingConsensusStructureMultiThreaded(double[][] basePairProb, double [] singleBaseProb)
    {
    	return performPosteriorDecodingMultiThreaded(basePairProb, singleBaseProb).pairedWith;
    }
    
    /**
     * Returns the posterior-decoding consensus structure.
     * @param basePairProb a NxN matrix of base-pairing probabilities. Assumes that sum(row) <= 1, in order to calculate the single base probabilities.
     * @return the posterior-decoding consensus structure.
     */
    public int [] getPosteriorDecodingConsensusStructureMultiThreaded(double[][] basePairProb)
    {
    	return performPosteriorDecodingMultiThreaded(basePairProb).pairedWith;
    }

    /**
     * A recursive method that fills the dynamic programming matrix for generating the posterior decoding structure.
     * @param basePairProb a NxN matrix of base-pairing probabilities.
     * @param singleBaseProb a vector of length N representing the probability that a base at a position is unpaired.
     * @param eMatrix the dynamic programming matrix for the posterior-decoding.
     * @param i the start position of the window.
     * @param j the end position of the window.
     * @param pairedWith an array of paired positions. Where (i, array[i]) represents a pairing between nucleotides (i+1, array[i]), if array[i] = 0, then (i+1) is unpaired.
     * @return the value of the eMatrix at position (i, j).
     */
    private static double recursePosteriorDecoding(double[][] basePairProb, double[] singleBaseProb, double[][] eMatrix, int i, int j, int[] pairedWith) {
        if (i > j) {
            return 0;
        }

        if (eMatrix[i][j] != RNAFoldingTools.emptyValue) {
            return eMatrix[i][j];
        }

        if (i == j) {
            eMatrix[i][j] = singleBaseProb[i];
            return eMatrix[i][j];
        }

        double u1 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, i + 1, j, pairedWith) + singleBaseProb[i];
        double p1 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, i + 1, j - 1, pairedWith) + basePairProb[i][j];
        double u2 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, i, j - 1, pairedWith) + singleBaseProb[j];
        double p2 = 0;
        for (int k = i; k < j; k++) {
            p2 = Math.max(p2, recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, i, k, pairedWith) + recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, k + 1, j, pairedWith));
        }

        eMatrix[i][j] = Math.max(u1, Math.max(p1, Math.max(u2, p2)));

        if (p1 > u1 && p1 > u2 && p1 > p2 && pairedWith != null) {
            pairedWith[i] = j + 1;
            pairedWith[j] = i + 1;
        }

        return eMatrix[i][j];
    }

    /**
     * Given a String array of dot-bracket structures, fills a vector which
     * counts the number of times a specific nucleotide position is unpaired.
     * For testing purposes only.
     * @param dotBracketStructures
     * @return an array of unpaired counts.
     */
    public static double[] getSingleBaseCount(String[] dotBracketStructures) {
        double[] singleBaseCount = new double[dotBracketStructures[0].length()];
        for (int i = 0; i < dotBracketStructures.length; i++) {
            String structure = dotBracketStructures[i];
            for (int j = 0; j < structure.length(); j++) {
                if (structure.charAt(j) == '.') // if unpaired
                {
                    singleBaseCount[j]++;
                }
            }
        }
        return singleBaseCount;
    }

    /**
     * Given an String array of dot-bracket structures, fills a matrix which 
     * counts the number of times a pair of nucleotides is base-paired in each structure.
     * For testing purposes only.
     * @param dotBracketStructures
     * @return a double matrix containing base-pairing counts.
     */
    public static double[][] getBasePairCountMatrix(String[] dotBracketStructures) {
        double[][] basePairMatrix = new double[dotBracketStructures[0].length()][dotBracketStructures[0].length()];

        for (int i = 0; i < dotBracketStructures.length; i++) {
            Stack<Integer> stack = new Stack<Integer>();
            String structure = dotBracketStructures[i];
            for (int j = 0; j < structure.length(); j++) {
                if (structure.charAt(j) == '(') {
                    stack.push(j);
                } else if (structure.charAt(j) == ')') {
                    int x = stack.pop();
                    int y = j;

                    basePairMatrix[x][y]++;
                    basePairMatrix[y][x]++;
                }
            }
        }

        return basePairMatrix;
    }

    /**
     * A helper method which prints double matrices.
     * For testing purposes only.
     * @param matrix 
     */
    public static void printMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(pad(matrix[i][j] + "", 4) + "  ");
            }
            System.out.println();
        }
    }

    /**
     * A helper method which pads or truncates a string to a specific length.
     * For testing purposes only.
     * @param s
     * @param length
     * @return the padded/truncated string.
     */
    public static String pad(String s, int length) {
        String ret = s;
        for (int i = s.length(); i < length; i++) {
            ret += " ";
        }
        return ret.substring(0, length);
    }

    public RNAFoldingTools() {
        //test();
    }

    public void test() {
        int repeats = 100;
        int size = 25;
        while (true) {
            System.out.println("-------------------------------------------------------------");
            size += 25;
            String s = "((";
            for (int i = 0; i < size; i++) {
                s += ".";
            }
            s += "))";
            String[] structures = {s};

            double[][] basePairCount = getBasePairCountMatrix(structures);
            double[] singleBaseCount = getSingleBaseCount(structures);

            //getPosteriorDecodingConsensusStructure(basePairCount, singleBaseCount);
            long startTime1 = System.currentTimeMillis();
            int divisionSize = 0;
            for (int i = 0; i < repeats; i++) {
                //getPosteriorDecodingConsensusStructure(basePairCount, singleBaseCount);
                MultiThreadedPosteriorDecoding m = new MultiThreadedPosteriorDecoding(basePairCount, singleBaseCount);
                m.compute();
                divisionSize = m.divisionSize;
            }
            long endTime1 = System.currentTimeMillis();
            long elapsed1 = endTime1 - startTime1;
            long startTime2 = System.currentTimeMillis();
            for (int i = 0; i < repeats; i++) {
                getPosteriorDecodingConsensusStructure(basePairCount, singleBaseCount);
            }
            long endTime2 = System.currentTimeMillis();
            long elapsed2 = endTime2 - startTime2;

            double rate1 = (double) repeats / (double) elapsed1;
            double rate2 = (double) repeats / (double) elapsed2;
            double ratio = rate1 / rate2;
            System.out.println((size + 4) + "\t" + elapsed1 + "\t" + elapsed2 + "\t" + divisionSize + "\t" + rate1 + "\t" + rate2 + "\t" + ratio);
            if (elapsed2 > 5000) {
                repeats = Math.max(1, repeats / 2);
            }
        }
    }

    /**
     * A class, which given a base-pairing probability matrix and an array 
     * representing the probabilities of a single nucleotides being unpaired, 
     * performs a multi-threaded posterior-decoding and returns the MPD consensus structure.
     */
    public class MultiThreadedPosteriorDecoding {

        double[][] basePairProb;
        double[] singleBaseProb;
        int[] pairedWith;
        final Integer lock = new Integer(0);
        final Integer waitLock = new Integer(0);
        int maxThreads = Runtime.getRuntime().availableProcessors();
        int threadsUsed = 0;
        int currentBlockY = 0;
        int currentX = 0;
        int currentY = 0;
        int length;
        int startingDivisionSize;
        int divisionSize;
        int blockIncrement;
        double[][] eMatrix;

        public MultiThreadedPosteriorDecoding(double[][] basePairCount, double[] singleBaseCount) {
            this.length = singleBaseCount.length;
            this.pairedWith = new int[this.length];
            this.basePairProb = basePairCount;
            this.singleBaseProb = singleBaseCount;
            //this.divisionSize = Math.max(length / 75, Math.min(length / maxThreads, maxThreads * maxThreads));
            this.divisionSize = Math.max(this.length / this.maxThreads, 10);
            this.blockIncrement = this.divisionSize;

            // this.divisionSize  = Math.max(this.length / this.maxThreads / 5,40);
            this.startingDivisionSize = this.divisionSize;
            //System.out.println(divisionSize);
            this.eMatrix = new double[length][length];
            for (int i = 0; i < eMatrix.length; i++) {
                for (int j = 0; j < eMatrix.length; j++) {
                    eMatrix[i][j] = RNAFoldingTools.emptyValue;
                }
            }
            compute();
        }

        /**
         * Initializes the posterior-decoding procedure by launching a specified
         * number of threads.
         * @see MultiThreadedPosteriorDecoding#computeNextSection() 
         * @return the value of the dynamic programming matrix at (0, N).
         */
        public double compute() {
            try {
                synchronized (waitLock) {

                    for (int i = 0; i < this.maxThreads; i++) {
                        computeNextSection();
                    }
                    waitLock.wait();
                }
            } catch (InterruptedException ex) {
                Logger.getLogger(RNAFoldingTools.class.getName()).log(Level.SEVERE, null, ex);
            }

            // recurse one last time *just* in case missed the 0,N case,
            return RNAFoldingTools.recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, 0, eMatrix.length - 1, pairedWith);
        }

        /**
         * Launches new threads until no more sections are available.
         * The thread class itself recursively calls computeNextSection() each time it completes,
         * ensuring that the number of threads running at any one time is less than number of
         * available processing cores.
         */
        public void computeNextSection() {
            synchronized (lock) {
                Pair p = getNextSection();
                if (p.x != -1) {
                    new PosteriorDecodingThread(p.x, p.y).start();
                }
            }
        }

        /**
         * Generates a sequence of coordinates of windows on which posterior-decoding can
         * be performed independently in parallel.
         * @return the positions of the next window to perform posterior-decoding.
         */
        public Pair getNextSection() {
            Pair p = getNextSectionRecursive();
            if (p.x == -1) {
                return p;
            } else {
                return new Pair(Math.min(p.x, this.length - 1), Math.min(p.y, this.length - 1));
            }
        }

        /**
         * @see MultiThreadedPosteriorDecoding#getNextSection() 
         */
        private Pair getNextSectionRecursive() {
            int newX = this.currentX;
            int newY = this.currentY;

            synchronized (lock) {
                if (newX == -1 && newY == -1) {
                    return new Pair(-1, -1);
                } else if (newX == 0 && newY == 0) {
                    newX = 0;
                    newY = divisionSize;
                } else {
                    newX += divisionSize;
                    newY += divisionSize;


                    if (currentBlockY >= this.length) {
                        newX = -1;
                        newY = -1;
                    } else if (newX >= this.length) {
                        currentBlockY += this.blockIncrement;

                        //int j = 
                        //
                        //divisionSize = (int)Math.max(Math.sqrt(2*Math.pow(length - currentBlockY,2))/ (double)maxThreads, 10);
                        divisionSize = (int) Math.max((double) (length - currentBlockY) / (double) maxThreads, 10);
                        //System.out.println(currentBlockY+"\t"+divisionSize);
                        //System.out.println(length+"\t"+currentBlockY+"\t"+divisionSize);
                        // start back at the first row of the matrix
                        //System.out.println(length+"\t"+startingDivisionSize);                        
                        //this.divisionSize = (int)Math.max((double) this.divisionSize - ((double)this.divisionSize / ((double)this.length / (double)this.divisionSize)), 1);
                        //System.out.println(this.divisionSize);

                        newX = 0;
                        newY = currentBlockY + divisionSize;
                    }
                }
            }

            this.currentX = newX;
            //this.currentX = this.currentY;
            /*
             * if(newX != 0) { this.currentX = this.currentY;
            }
             */
            this.currentY = newY;

            if (newY - divisionSize > length) {
                return getNextSection();
            }

            return new Pair(newX, newY);
        }

        /**
         * A single thread which performs posterior-decoding on a specified window.
         */
        class PosteriorDecodingThread extends Thread {

            int x;
            int y;

            public PosteriorDecodingThread(int x, int y) {
                this.x = x;
                this.y = y;
            }

            public void run() {
                RNAFoldingTools.recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, this.x, this.y, pairedWith);
                computeNextSection();

                if (this.x == 0 && this.y == MultiThreadedPosteriorDecoding.this.length - 1) {
                    // complete, so notify
                    synchronized (waitLock) {
                        waitLock.notify();
                    }
                }
            }
        }
    }
    
    /**
     * Returns a dot bracket string representation given an array paired sites.
     * @param pairedSites
     * @return a dot bracket string representation given an array paired sites.
     */
    public static String getDotBracketStringFromPairedSites(int [] pairedSites)
    {
    	String dbs = "";
    	for(int i = 0 ; i < pairedSites.length ; i++)
    	{
    		if(pairedSites[i] == 0)
    		{
    			dbs += ".";
    		}
    		else
    		if(pairedSites[i] > i+1)
    		{
    			dbs += "(";
    		}
    		else
    		{
    			dbs += ")";
    		}
    	}
    	return dbs;
    }
    
    /**
     * Given a dot bracket string representation, returns an array of paired sites.
     * @param dotBracketStructure a dot bracket string representation of a structure.
     * @return an array of paired sites.
     * @see #getDotBracketStringFromPairedSites(int[])
     */
    public static int [] getPairedSitesFromDotBracketString(String dotBracketStructure)
    {
    	int [] pairedSites = new int[dotBracketStructure.length()];
        Stack<Integer> stack = new Stack<Integer>();
        for (int j = 0; j < dotBracketStructure.length(); j++) {
            if (dotBracketStructure.charAt(j) == '(') {
                stack.push(j);
            } else if (dotBracketStructure.charAt(j) == ')') {
                int x = stack.pop();
                int y = j;
                pairedSites[x] = y + 1;
                pairedSites[y] = x + 1;
            }
        }
        return pairedSites;
    }
    
    /**
     * Given a dot bracket string representation, returns an array of paired sites.
     * @param dotBracketStructure a dot bracket string representation of a structure.
     * @param openBracket the opening bracket.
     * @param closeBracket the closing bracket.
     * @return an array of paired sites.
     * @see #getDotBracketStringFromPairedSites(int[])
     */
    public static int [] getPairedSitesFromDotBracketString(String dotBracketStructure, char openBracket, char closeBracket)
    {
    	int [] pairedSites = new int[dotBracketStructure.length()];
        Stack<Integer> stack = new Stack<Integer>();
        for (int j = 0; j < dotBracketStructure.length(); j++) {
            if (dotBracketStructure.charAt(j) == openBracket) {
                stack.push(j);
            } else if (dotBracketStructure.charAt(j) == closeBracket) {
                int x = stack.pop();
                int y = j;
                pairedSites[x] = y + 1;
                pairedSites[y] = x + 1;
            }
        }
        return pairedSites;
    }
}
