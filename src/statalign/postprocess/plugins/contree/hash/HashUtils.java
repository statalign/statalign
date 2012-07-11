package statalign.postprocess.plugins.contree.hash;

import java.util.Random;

/**
 * Details about this class can be found in papers on consensus tree construction algorithms:
 * - An Experimental Analysis of Consensus Tree Algorithm [..]. Seung-Jin Sul et al. 2009.
 * - A Linear-time Majority Tree Algorithm. Amenta et al. 2004.
 * 
 * @author eiriksson
 */
public class HashUtils {

    // Constants

    private static Random random = new Random();

    // Variables

    public int m1, m2;
    public int[] a1, a2;

    // Functions

    /** Initializes the hash utilities. */
    public void initialize(int noOfTaxa, int noOfTrees, int c, long seed) {        
        // Chooses the size of hash table and the
        m1 = getNextPrime(noOfTaxa * noOfTrees);
        m2 = getNextPrime(noOfTaxa * noOfTrees * c);

        // Sets the seed.
        random.setSeed(seed);

        // Initializes the arrays which the universal hash functions uses.
        a1 = new int[noOfTaxa];
        a2 = new int[noOfTaxa];
        for (int i = 0; i < noOfTaxa; i++) {
            a1[i] = random.nextInt(m1);
            a2[i] = random.nextInt(m2);
        }
    }

    // Static functions

    /** Determines the next prime above or equal to n. */
    public static int getNextPrime(int n) {
        if ((n & 1) == 0) {
            n++;
        }
        while (!isPrime(n)) n += 2;
        return n;
    }

    /** Determines if a number is a prime or not. */
    public static boolean isPrime(long n) {
        boolean prime = true;
        for (long i = 3; i <= Math.sqrt(n); i += 2)
            if (n % i == 0) {
                prime = false;
                break;
            }
        if (((n & 1) == 1 && prime && n > 2) || n == 2) {
            return true;
        } else {
            return false;
        }
    }

    /** For debugging purposes only. */
    public static void main(String[] args) {
        assert getNextPrime(8) == 11;
        assert getNextPrime(9) == 11;
        assert getNextPrime(11) == 11;
        assert getNextPrime(75884) == 75913;
        assert getNextPrime(9124961) == 9124967;
    }

}
