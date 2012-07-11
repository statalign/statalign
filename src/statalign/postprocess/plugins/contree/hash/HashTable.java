package statalign.postprocess.plugins.contree.hash;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.LinkedList;

public class HashTable {

    // Variables

    /** The main hash table. */
    private ArrayList<HashEntry>[] hashTable;

    // Functions

    @SuppressWarnings({"unchecked"})
    public HashTable(int size) {
        hashTable = (ArrayList<HashEntry>[]) new ArrayList[size]; // Jeez - Java generics => horrible.
        for (int i = 0; i < size; i++) {
            hashTable[i] = new ArrayList<HashEntry>();
        }
    }

    /**
     * Puts a bi-partition into the hash table. If this bi-partition becomes a majority partition
     * this function will also put that into the <code>partitions</code> linked list.
     * @param partition the bi-partition.
     * @param tableKey the key that indexes the actual hash table.
     * @param bucketKey the bucket key which is kept in each entry.
     * @param resRate how many trees a bi-partition has to be in, to be a majority bi-partition.
     * @param partitions a {@link LinkedList} containing the majority partitions.
     */
    public void put(BitSet partition, double edgeLength,
                    int tableKey, int bucketKey, double resRate,
                    LinkedList<HashEntry> partitions) {
        if (hashTable[tableKey].isEmpty()) {                        // No entry exists for this key (1): Create.
            HashEntry entry = new HashEntry(bucketKey, partition, edgeLength);
            if (1.0 > resRate && !entry.isMajority) {   // A new partition of interest!
                entry.isMajority = true;  // is majority if is of interest at this point...!
                partitions.add(entry);
            }
            hashTable[tableKey].add(entry);
        } else {
            boolean found = false;
            for (HashEntry entry : hashTable[tableKey]) {           // Searches the bucket.
                if (entry.hashKey2 == bucketKey) {                  // An entry was found.
                    entry.count++;                                      // Increase the occurrence of the entry.
                    entry.edgeLengthsSum += edgeLength;
                    if ((double)entry.count > resRate && !entry.isMajority) {   // A new partition of interest!
                        entry.isMajority = true;  // is majority if is of interest at this point...!
                        partitions.add(entry);
                    }
                    found = true;
                    break;
                }
            }
            if (!found) {                                               // No entry exists for this key (2): Create.
                HashEntry entry = new HashEntry(bucketKey, partition, edgeLength);
                if (1.0 > resRate && !entry.isMajority) {   // A new partition of interest!
                    entry.isMajority = true;  // is majority if is of interest at this point...!
                    partitions.add(entry);
                }
                hashTable[tableKey].add(entry);
            }
        }
    }

}
