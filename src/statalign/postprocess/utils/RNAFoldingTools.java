package statalign.postprocess.utils;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Michael
 */
public class RNAFoldingTools {
	
    public static void main(String[] args) {
    	ArrayList<String> sequences = new ArrayList<String>();
    	sequences.add("ACT--CC-");
    	sequences.add("ACTC-CCG");
    	sequences.add("ACTCTCCG");
    	System.out.println(RNAFoldingTools.getReferenceSequence(sequences, 9));
    	
    	
    	// Generate an example "Base pairing probability" matrix
    	String [] structures = {"(((..(..))))(..)", 
    							"(((..(..))))(..)", 
    							"(((..(..))))(..)"};
    	double[][] basePairCount = getBasePairCountMatrix(structures);
    	double[] singleBaseCount = getSingleBaseCount(structures);
    	
    	// Perform the posterior decoding
    	RNAFoldingTools rnaTools = new RNAFoldingTools();
    	//int [] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(basePairCount); // fails for this example, because not base-pairing probability matrix
    	//int [] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(basePairCount, singleBaseCount); // this works
    	
    	int [] pairedSites = RNAFoldingTools.getPosteriorDecodingConsensusStructure(basePairCount, singleBaseCount);
    	
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
    
    public static double [] getSingleBaseProb(double[][] basePairProb)
    {
    	double [] singleBaseProb = new double[basePairProb.length];
    	for(int i = 0 ; i < basePairProb.length ; i++)
    	{
    		singleBaseProb[i] = 1;
    		for(int j = 0 ; j < basePairProb[0].length ; j++)
        	{
    			singleBaseProb[i] -= basePairProb[i][j];
        	}    		
    	}
    	return singleBaseProb;
    }
    
    public static int [] getPosteriorDecodingConsensusStructure(float[][] basePairProb) {    
    	return getPosteriorDecodingConsensusStructure(getDoubleMatrix(basePairProb));
    }
    
    public static int [] getPosteriorDecodingConsensusStructure(double[][] basePairProb) {
    	double [] singleBaseProb = new double[basePairProb.length];
    	for(int i = 0 ; i < basePairProb.length ; i++)
    	{
    		singleBaseProb[i] = 1;
    		for(int j = 0 ; j < basePairProb[0].length ; j++)
        	{
    			singleBaseProb[i] -= basePairProb[i][j];
        	}    		
    	}
    	return getPosteriorDecodingConsensusStructure(basePairProb, singleBaseProb);
    }

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
        int[][] S = new int[basePairProb.length][basePairProb[0].length];
        recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, 0, eMatrix.length-1);

        //printMatrix(eMatrix);
        //System.out.println();
        //printMatrix(S);
        writeMatrix(S, new File("e.matrix"));
        writeMatrix(S, new File("s.matrix"));
        traceBack(S, 0, eMatrix.length - 1, pairedWith);
        
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

    	/*for(int i = 0 ; i < singleBaseProb.length ; i++)
    	{
    		System.out.println((i+1)+"\t"+singleBaseProb[i]);
    	}*/
    	
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
    
    
    public int [] getPosteriorDecodingConsensusStructureMultiThreaded(float[][] basePairProb) {    
    	return getPosteriorDecodingConsensusStructureMultiThreaded(getDoubleMatrix(basePairProb));
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
        double p1 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, i + 1, j - 1, pairedWith) + basePairProb[i][j]; // * 2
        double u2 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, i, j - 1, pairedWith) + singleBaseProb[j];
        double p2 = 0;
        for (int k = i; k < j; k++) {
        	// remember K
            p2 = Math.max(p2, recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, i, k, pairedWith) + recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, k + 1, j, pairedWith));
        }

        eMatrix[i][j] = Math.max(u1, Math.max(p1, Math.max(u2, p2)));

        if (p1 > u1 && p1 > u2 && p1 > p2 && pairedWith != null) {
        	
        	// if(pairedWith[i] == 0 && pairedWith[j] == 0)
        	 //{
	            pairedWith[i] = j + 1;
	            pairedWith[j] = i + 1;
        	 //}
	            System.out.println("B"+i+"\t"+j);
        }

        return eMatrix[i][j];
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
    private static double recursePosteriorDecoding(double[][] basePairProb, double[] singleBaseProb, double[][] eMatrix, int [][] S, int i, int j) {
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

        
        double u1 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, i + 1, j) + singleBaseProb[i];
        double p1 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, i + 1, j - 1) + 2*basePairProb[i][j]; // * 2
        double u2 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, i, j - 1) + singleBaseProb[j];
        double p2 = 0;
        int max_k = i;
        for (int k = i; k < j; k++) {
        	double p2k =  recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, i, k) + recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, k + 1, j);
        	// remember K
            if(p2k > p2)
            {
            	p2 = p2k;
            	max_k = k;
            }
        }

        /*
        double max = 0;
        if(p1 > u1 && p1 > p2 && p1 > u2)
        {
        	max = p1;
        	S[i][j] = -3;
        }
        else
        if(u1 > p2 && u1 > u2)
        {
        	max = u1;
        	S[i][j] = -1;
        }
        else
        if(p2 > u2)
        {
        	max = p2;
        	S[i][j] = max_k;
        }
        else
        {
        	max = u2;
        	S[i][j] = -2;
        }*/
        
       
        double max = 0;
        if(u1 > p1 && u1 > p2 && u1 > u2)
        {
        	max = u1;
        	S[i][j] = -1;
        }
        else
    	if(u2 > p1 && u2 > p2)
        {
        	max = u2;
        	S[i][j] = -2;
        }
        else
        if(p1 > p2)
        {
        	max = p1;
        	S[i][j] = -3; // paired
        }
        else
        {
        	max = p2;
        	S[i][j] = max_k;
        }
        
        eMatrix[i][j] = max;
        
        return eMatrix[i][j];
    }
    
    public static void traceBack(int [][] S, int i, int j, int [] pairedWith)
    {
    	if(i >= j)
    	{
    		// do nothing
    	}
    	else
    	if(S[i][j] == -3)
    	{
    		pairedWith[i] = j + 1;
    		pairedWith[j] = i + 1;
    		traceBack(S, i+1, j-1, pairedWith);    
    	}
    	else
    	{
    		traceBack(S, i, S[i][j], pairedWith);
    		traceBack(S, S[i][j]+1, j, pairedWith);
    	}
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
     * A helper method which writes double matrices to files.
     * For testing purposes only.
     * @param matrix 
     */
    public static void writeMatrix(double[][] matrix, File file) {
    	DecimalFormat df = new DecimalFormat("0.0000E0");
    	try
    	{
	    	BufferedWriter buffer = new BufferedWriter(new FileWriter(file));
	        for (int i = 0; i < matrix.length; i++) {
	            for (int j = 0; j < matrix[0].length; j++) {
	            	buffer.write(pad(df.format(matrix[i][j]) + "", 10) + "  ");
	            }
	            buffer.newLine();
	        }
	        buffer.close();
    	}
    	catch(IOException ex)
    	{
    		ex.printStackTrace();
    	}
    }
    
    /**
     * A helper method which writes integer matrices to files.
     * For testing purposes only.
     * @param matrix 
     */
    public static void writeMatrix(int[][] matrix, File file) {
    	try
    	{
	    	BufferedWriter buffer = new BufferedWriter(new FileWriter(file));
	        for (int i = 0; i < matrix.length; i++) {
	            for (int j = 0; j < matrix[0].length; j++) {
	            	buffer.write(pad(matrix[i][j] + "", 4) + "  ");
	            }
	            buffer.newLine();
	        }
	        buffer.close();
    	}
    	catch(IOException ex)
    	{
    		ex.printStackTrace();
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
        int[][] S;

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
            this.S = new int[length][length];
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
            RNAFoldingTools.recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, 0, eMatrix.length - 1, pairedWith);
            
            RNAFoldingTools.traceBack(S, 0, eMatrix.length-1, pairedWith);
            return eMatrix[0][eMatrix.length-1];
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
                RNAFoldingTools.recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, this.x, this.y);
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
            	if(!stack.empty())
            	{
	                int x = stack.pop();
	                int y = j;
	                pairedSites[x] = y + 1;
	                pairedSites[y] = x + 1;
            	}
            }
        }
        return pairedSites;
    }
    
    public static String getReferenceSequence(ArrayList<String> sequences, int refLength)
    {
    	int length = 0;
    	for(int i = 0 ; i < sequences.size() ; i++)
    	{
    		length = Math.max(length, sequences.get(i).length());
    	}
    	
    	int [] counts = new int[length];
    	int maxCount = 0;
    	boolean [] countUsed = new boolean[length];
    	for(int i = 0 ; i < sequences.size() ; i++)
    	{
    		String seq = sequences.get(i);
    		for(int j = 0 ; j < seq.length() ; j++)
    		{
    			if(seq.charAt(j) != '-')
    			{
    				counts[j] += 1;
    				maxCount = Math.max(maxCount, counts[j]);
    			}
    		}
    	}
    	
    	String referenceString = "";
    	int n = 0;
    	for(int k = maxCount ; k >= 0 ; k--)
    	{
	    	for(int i = 0 ; i < counts.length ; i++)
	    	{
	    		if(!countUsed[i] && counts[i] == k)
	    		{
	    			countUsed[i] = true;
	    			n++;
	    		}
	    		
	    		if(n == refLength)
	    		{
	    			break;
	    		}
	    	}
	    	
	    	if(n == refLength)
	    	{
	    		break;
	    	}
    	}
    	
    	String ref = "";
    	for(int i = 0 ; i < countUsed.length ; i++)
    	{
    		if(countUsed[i])
    		{
    			ref += "N";
    		}
    		else
    		{
    			ref += "-";
    		}
    	}
    	
    	while(ref.length() < refLength)
    	{
    		ref += "X";
    	}
    	
    	return ref;
    	
    }
    
    public static String getDotBracketStringFromCtFile(File ctFile)
    {
    	return RNAFoldingTools.getDotBracketStringFromPairedSites(RNAFoldingTools.getPairedSitesFromCtFile(ctFile));
    }
    
    public static String getSequenceByName(String seqName, ArrayList<String> sequences, ArrayList<String> sequenceNames)
    {
    	for(int i = 0 ; i < sequences.size() ; i++)
    	{
    		if(sequenceNames.get(i).equals(seqName))
    		{
    			return sequences.get(i);
    		}
    	}
    	System.err.println("Could not find ref seq: "+seqName);
    	return sequences.get(0);
    }
    
    public static int [] getPairedSitesFromCtFile(File ctFile)
    {
    	try
    	{
	    	 BufferedReader buffer = new BufferedReader(new FileReader(ctFile));
	
	         String textline = null;
	
	         int [] pairedSites = null;
	         while ((textline = buffer.readLine()) != null) {
	             String[] split = textline.trim().split("(\\s)+");
	             int length = Integer.parseInt(split[0]);
	             pairedSites = new int[length];
	             //String sequence = "";
	             for (int i = 0; i < length && (textline = buffer.readLine()) != null; i++) {
	                 String[] split2 = textline.trim().split("(\\s)+");
	                 pairedSites[Integer.parseInt(split2[0])-1] =  Integer.parseInt(split2[4]);
	             }
	         }
	
	         buffer.close();
	         
	         return pairedSites;
    	}
    	catch(IOException ex)
    	{
    		ex.printStackTrace();
    	}
    	
    	return null;
    }
    
    public static int [] getPairedSitesFromDBNStringFile(File dbnFile)
    {
    	try
    	{
	    	 BufferedReader buffer = new BufferedReader(new FileReader(dbnFile));
	    	 String textline = buffer.readLine();
	    	 buffer.close();
	    	 System.out.println(textline);
	    	 return getPairedSitesFromDotBracketString(textline);
    	}
    	catch(IOException ex)
    	{
    		ex.printStackTrace();
    	}
    	return null;
    }
    
    public static double [][] loadMatrix(File bpFile)
    {
    	//File bpFile = new File("C:/Users/Michael/Dropbox/RNA and StatAlign/TestRNAData/TestRNAData1.dat.bp");		
		double [][] bpMatrix = null;
		try
		{
			BufferedReader buffer = new BufferedReader(new FileReader(bpFile));
			String textline = buffer.readLine();
			int length = textline.split("(\\s)+").length;
			bpMatrix = new double[length][length];

			for(int i = 0 ; i < length ; i++)
			{
				String [] split = textline.split("(\\s)+");
				for(int j = 0 ; j < length ; j++)
				{
					bpMatrix[i][j] = Double.parseDouble(split[j]);
				}
				textline = buffer.readLine();
			}
			
			buffer.close();

		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
		
		return bpMatrix;
		
    }
    
    public static double [][] getDoubleMatrix(float [][] matrix)
    {
    	double [][] doubleMatrix = new double[matrix.length][matrix[0].length];
    	for(int i = 0 ; i < matrix.length ; i++)
    	{
    		for(int j = 0 ; j < matrix[0].length ; j++)
    		{
    			doubleMatrix[i][j] = matrix[i][j];
    		}
    	}
    	return doubleMatrix;
    }
    
    public static float [][] getFloatMatrix(double [][] matrix)
    {
    	float [][] floatMatrix = new float[matrix.length][matrix[0].length];
    	for(int i = 0 ; i < matrix.length ; i++)
    	{
    		for(int j = 0 ; j < matrix[0].length ; j++)
    		{
    			floatMatrix[i][j] = (float)matrix[i][j];
    		}
    	}
    	return floatMatrix;
    }
    
    public static void loadFastaSequences(File file, ArrayList<String> sequences, ArrayList<String> sequenceNames) {
        loadFastaSequences(file, sequences, sequenceNames, Integer.MAX_VALUE);
    }

    public static void loadFastaSequences(File file, ArrayList<String> sequences, ArrayList<String> sequenceNames, int max) {
        try {
            BufferedReader buffer = new BufferedReader(new FileReader(file));
            String textline = null;
            String sequence = "";
            int n = 0;
            boolean maxReached = false;
            while ((textline = buffer.readLine()) != null) {
                if(maxReached && textline.startsWith(">"))
                {
                    break;
                }
                
                if (textline.startsWith(">")) {
                    n++;
                    if (n >= max) {
                        maxReached = true;
                    }

                    sequenceNames.add(textline.substring(1));
                    if (!sequence.equals("")) {
                        sequences.add(sequence.toUpperCase());
                        sequence = "";
                    }
                } else {
                    sequence += textline.trim();
                }

            }
            buffer.close();
            if (!sequence.equals("")) {
                sequences.add(sequence);
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
	
	public static double calculatePPfoldReliabilityScore(int [] pairedSites, double [][] basePairProb)
	{
		double [] singleBaseProb = RNAFoldingTools.getSingleBaseProb(basePairProb);
		double ppfoldReliablityScore = 0;
		for(int i = 0 ; i < pairedSites.length ; i++)
		{
			if(pairedSites[i] == 0)
			{
				ppfoldReliablityScore += singleBaseProb[i];
			}
			else
			{
				ppfoldReliablityScore += basePairProb[i][pairedSites[i]-1];
			}
		}

		return ppfoldReliablityScore/((double)pairedSites.length);
	}
	
	public static double calculatePairsOnlyReliabilityScore(int [] pairedSites, double [][] basePairProb)
	{
		double ppfoldReliablityScore = 0;
		double pairs = 0;
		for(int i = 0 ; i < pairedSites.length ; i++)
		{
			if(pairedSites[i] != 0)
			{
				ppfoldReliablityScore += basePairProb[i][pairedSites[i]-1];
				pairs++;
			}
		}
		
		if(pairs == 0)
		{
			return 1;
		}

		return ppfoldReliablityScore/pairs;
	}
	
	public static double calculatePairsOnlyReliabilityScore(int [] pairedSites, double [][] basePairProb, ArrayList<Double> weights)
	{
		double ppfoldReliablityScore = 0;
		double pairs = 0;
		double weightSum = 0;
		for(int i = 0 ; i < pairedSites.length ; i++)
		{
			if(pairedSites[i] != 0)
			{
				double weighting = weights.get(i)*weights.get(pairedSites[i]-1);
				ppfoldReliablityScore += basePairProb[i][pairedSites[i]-1]*weighting;
				weightSum += weighting;
				pairs++;
			}
		}
		
		if(pairs == 0)
		{
			return 1;
		}

		return ppfoldReliablityScore/pairs;
	}
	
	public static void writeToFile(File f, String s, boolean append)
	{
		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(f, append));
			buffer.write(s+"\n");
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	public static void saveCtFile(File outFile, int [] pairedSites, String header, String sequence)
	{
		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile));
			buffer.write(header+"\n");
			for(int i = 0 ; i < pairedSites.length ; i++)
			{
				buffer.write((i+1)+sequence.charAt(i)+"\t"+i+"\t"+(i+2)+"\t"+pairedSites[i]+"\t"+(i+1)+"\n");
			}
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	public static void saveDotBracketFile(File outFile, int [] pairedSites, String header, String sequence)
	{
		try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter(outFile));
			buffer.write(">"+header+"\n");
			buffer.write(sequence+"\n");
			buffer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites)+"\n");
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	public static boolean isRNAalignment(ArrayList<String> sequences)
	{
		double countRNA = 0;
		double countNonRNA = 0;
		for(int i = 0 ; i < sequences.size() ; i++)
		{
			String s = sequences.get(i).toUpperCase();
			for(int j = 0 ; j < s.length() ; j++)
			{
				switch(s.charAt(j))
				{
				case 'A':
					countRNA++;
					break;
				case 'C':
					countRNA++;
					break;
				case 'G':
					countRNA++;
					break;
				case 'T':
					countRNA++;
					break;
				case 'U':
					countRNA++;
					break;
				case '-':
					break;
				default:
					countNonRNA++;
				}
			}			
		}
		
		double ratio = countRNA / countNonRNA;		
		return ratio > 0.5;
	}
}
