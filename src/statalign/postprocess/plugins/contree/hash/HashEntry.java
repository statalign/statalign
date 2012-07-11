package statalign.postprocess.plugins.contree.hash;

import java.util.BitSet;

public class HashEntry {

    /** The bi-partition relevant to this entry. */
    public BitSet partition;

    /** The bucket hash key relevant to this entry. */
    public int hashKey2;

    /** The number of times this bi-partition has been seen. */
    public int count;

    /** If this is a majority partition this will be set to true, otherwise false. */
    public boolean isMajority;

    /** The sum of the edge lengths */
    public double edgeLengthsSum;

    public double[] leaves;

    public HashEntry(int bucketKey) {
        this.hashKey2 = bucketKey;
        this.count = 1;
        isMajority = false;
    }

    public HashEntry(int bucketKey, BitSet partition, double edgeLength) {
        this(bucketKey);
        this.partition = partition;
        this.edgeLengthsSum += edgeLength;
    }

}
