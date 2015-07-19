package statalign.postprocess.plugins.contree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import statalign.postprocess.plugins.TreeNode;
import statalign.postprocess.plugins.contree.hash.HashEntry;
import statalign.postprocess.plugins.contree.hash.HashTable;
import statalign.postprocess.plugins.contree.hash.HashUtils;
/**
 * The main thread for calculating consensus trees and networks.
 * Note that networks is dependent on the initial tree...
 * 
 * @author wood, eiriksson
 *
 */
public class CTMain {

    // Constants

    /** The default double collision constant - the larger, the better. */
    private final static int C = 1000;

    /** The default resolution rate in percentage. 
     * 33.34 - splits must occur in over a third of the trees for it to ALLWAYS be possible to draw in 2D, no matter what the input TREES are.
     * */
    private final static double DEFAULT_RES_PERCENTAGE = 33.34; 


    // Variables
    private int noOfSamples;                // Number of samples.
    private int noOfTaxa;                   // Number of taxa.
    private double resPercentage;           // The resolution rate percentage (0 < x <= 100). 
    private long seed;                      // The seed used in the random number generator.
    private TaxaMap taxa;                   // HashMap containing mapping from taxon names -> indices

    private HashUtils hashUtils;            // Holds stuff concerning the hashing.
    private HashTable hashTable;            // The hash table.

    private int noOfTrees;                  // Current number of trees.
    private int majorityThreshold;          // Declares in how many trees partitions need to appear to be a majority partition.
    private double interestThreshold;          // Declares in how many trees partitions need to appear to be of interest. 

    private LinkedList<HashEntry> partitions;// Holds the partitions of interest. 
    private int noOfPartitions;             // Holds the current number of partitions.

    private double[] leafEdgeLengths;		// Array of leaf edge lengths
    
    private TreeNode root; 	// The tree root of the input tree.

    public boolean initialized = false;
    // Functions
    /**
     * Initialise the class
     */
    public CTMain() {
        // Initializes some defaults (can be changed via setters).
        Random random = new Random();
        this.resPercentage = DEFAULT_RES_PERCENTAGE;
        this.seed = random.nextLong();

        // Object initialisation 
        partitions = new LinkedList<HashEntry>();
    }
    /**
     * Print some configuration information
     *
     */
    public void printConfig() {
        System.out.printf("Number of taxa: %d - ", noOfTaxa);
        System.out.printf("Number of trees: %d - ", noOfTrees);
        System.out.printf("Resolution: %d\n", interestThreshold); 
    }
    /**
     * Initialise a CTMain before usage with a first tree, going through and setting up the hash table etc.
     * 
     * @param root Root of the initial tree to initialise with
     * @param noOfSamples Current number of samples taken
     *
     */
    public void initialize(TreeNode root, int noOfSamples) {
        // Parameter initialisation
        this.noOfSamples = noOfSamples;
        noOfTrees = 0;
        // TaxaMap initialisation
        List<TreeNode> leaves = root.getLeaves();
        noOfTaxa = leaves.size();
        // Hash initialisation
        hashUtils = new HashUtils();
        hashUtils.initialize(noOfTaxa, noOfSamples, C, seed);
        hashTable = new HashTable(hashUtils.m1);
        //Taxamap initialisation
        taxa = new TaxaMap(noOfTaxa);
        for (int i = 0; i < leaves.size(); i++) {
            taxa.put(leaves.get(i).name, i);
        }
        leafEdgeLengths = new double[noOfTaxa];
        // Adds a single star partition, once and for all.
        BitSet star = new BitSet(noOfTaxa);
        star.flip(0, noOfTaxa);
        HashEntry entry = new HashEntry(-1, star, 0.0d);
        entry.count = noOfSamples + 1;
        partitions.add(entry);
        // Majority threshold initialisation
        updateInterestThreshold(); 
        initialized = true;
    } 
    /**
     * Simply update our interest threshold
     */
    // Update the Interest Threshold
    public void updateInterestThreshold() { 
    	interestThreshold = (double) ((double) noOfTrees * (resPercentage / 100.0d)); 
	}
    /**
     * Create partitions from an input tree - recursive so will be called by many nodes, beginning with root but with calculations actually beginning on leaves.
     *  
     * @param node Node in tree with above split to add to the hash table
     *
     */
	private BitSet createPartitions(TreeNode node) { 
		if (node.isLeaf()) { // Leaf node. 
        	int index = taxa.get(node.name);
        	// leaf simply has it's keys stored in the hash utilities
            node.addProperty("tableHashKey", hashUtils.a1[index]);
            node.addProperty("bucketHashKey", hashUtils.a2[index]);
            // Updates the edge length array.
            leafEdgeLengths[index] += node.edgeLength; 
            // Create a new partition from scratch to represent it
            BitSet partition = new BitSet(noOfTaxa); 
            partition.set(index);
            assert partition.cardinality() == 1 : "There should be exactly a single bit set."; 
            return partition;
		} else { // An internal node: Traverses the tree in post order. 
        	BitSet partition = new BitSet(noOfTaxa);
        	// Get the node's partition representation from its children
        	List<TreeNode> children = node.children;
        	for (TreeNode child : children) {
                partition.or(createPartitions(child)); 
        	}
        	// if this node is NOT to the right side of the root then add it... (Avoid adding splits twice for CNetworks)
        	if (node.parent != root || root.getRight() != node) { 
                noOfPartitions++; 
                long tableKey = 0; 
                long bucketKey = 0; 
                long tableKey2 = 0; 
                long bucketKey2 = 0; 
                // Calculate the hash keys for this partition
                for (TreeNode child : children) { 
                    tableKey += child.getIntProperty("tableHashKey");
                    bucketKey += child.getIntProperty("bucketHashKey");
                } 
                // if first is one then we need to store the flipped version so we store each split in one representation only.
                if (partition.get(0) == true) { 
                    // copy to a new partition that is the flipped version
                    BitSet partitionF = new BitSet(noOfTaxa); 
                    for(int l=0;l<noOfTaxa;l++){ 
                    	if(partition.get(l)==false)partitionF.set(l); 
                    } 
                    // calculate the hash keys for the flipped partition...
                    for (int k=0;k<noOfTaxa;k++) { 
                        if(partitionF.get(k)==true){ 
                        	tableKey2 += hashUtils.a1[k]; 
                        	bucketKey2 += hashUtils.a2[k]; 
                        } 
                    } 
                    // store the properties in the node
                    node.addProperty("tableHashKey", (int) (tableKey2 % hashUtils.m1));  
                    node.addProperty("bucketHashKey", (int) (bucketKey2 % hashUtils.m2));
                    if (noOfPartitions < noOfTaxa - 2) { // Avoids the addition of the star partition 
                            hashTable.put(partitionF, node.edgeLength,node.getIntProperty("tableHashKey"), node.getIntProperty("bucketHashKey"), interestThreshold, partitions); 
                    } 
                    // remember to still return the original partition with its appropriate keys...
                    node.addProperty("tableHashKey", (int) (tableKey % hashUtils.m1));  
                    node.addProperty("bucketHashKey", (int) (bucketKey % hashUtils.m2)); 
                    return partition; 
                } 
                // if first is zero then simply add the partition, recursively calculating the hash
                else{ 
                    node.addProperty("tableHashKey", (int) (tableKey % hashUtils.m1));  
                    node.addProperty("bucketHashKey", (int) (bucketKey % hashUtils.m2)); 
                    if (noOfPartitions < noOfTaxa - 2) { // Avoids the addition of the star partition 
                    	hashTable.put(partition, node.edgeLength, node.getIntProperty("tableHashKey"),node.getIntProperty("bucketHashKey"), interestThreshold, partitions); 
                    } 
                    return partition; 
                } 
        	} 
        	// still return the partition even if node was not added as was right of root...
        	return partition; 
		} 
	} 
	 /**
     * Add a new tree to the hash table
     * 
     * @param root Root of new tree to add
     *
     */
    public void addNewTree(TreeNode root) {
        // Updates the number of trees and the threshold.
        noOfTrees++;
        updateInterestThreshold();

        // Hashes the partitions of the trees.
        noOfPartitions = 0;
        this.root = root; 
        createPartitions(root);
    }
    /**
     * Create partitions from an input tree - recursive so will be called by many nodes, beginning with root but with calculations actually beginning on leaves.
     *  
     * @param partitions Partitions in the form of entries in the hash table
     * @param curInterestPercentage Current percentage of interest that we want splits to occur above to view in the network later
     *
     */
	private ArrayList<Cluster> constructClusters(LinkedList<HashEntry> partitions, double curInterestPercentage) {
		ArrayList<Cluster> clusters = new ArrayList<Cluster>();
		// Thresholds for the current run of the cluster builder (c.f. the threshold for the partitions list).
		double curInterestThreshold = (double) (noOfTrees * (curInterestPercentage / 100.0d));
		int majInterestThreshold = (int) (noOfTrees * (50.0 / 100.0d));
		for (Iterator<HashEntry> it = partitions.iterator(); it.hasNext();) {
			HashEntry entry = it.next();
			// Checks if this partition is still above threshold... partition.
			// If not: remove it - O(1).
			if ((double)entry.count <= curInterestThreshold) {
				// correct the isMajority flag to now refer to majority not if of interest or not..
				if (entry.count <= majInterestThreshold) {
					entry.isMajority = false;
				}
				it.remove();
				continue;
			}
			// Constructs clusters (list of TreeNode's for each set bit) from each partition.
			if ((double)entry.count > curInterestThreshold) {
				Cluster cluster = new Cluster();
				if (entry.count > majInterestThreshold) {
					cluster.isMajority = true;
				}
				cluster.aboveSplit = entry.partition;
				cluster.noOfOccurrences = entry.count;
				cluster.edgeLength = entry.edgeLengthsSum / entry.count;
				for (int i = 0; i < entry.partition.size(); i++) {
					if (entry.partition.get(i)) {
						TreeNode node = new TreeNode(taxa.getName(i));
						node.edgeLength = leafEdgeLengths[i] / noOfTrees;
						cluster.add(node);
					}
				}
				clusters.add(cluster);
			}
		}
		// Sort by number of taxa.
		// TODO: This might obviously be optimized a bit, e.g. with a PriorityQueue. - Eiriksson
		Collections.sort(clusters);
		return clusters;
	}
    /**
     * Constructs the majorityTree.  Used to make Network so remember to call this first!  
     */
    public CTree constructMajorityTree() {
        int numOfNode = 0;
        CTree tree = new CTree();
        // Creates the clusters.
        tree.clusters = constructClusters(partitions,resPercentage);


        // Begins by constructing the star tree.

        for (TreeNode node : tree.clusters.get(0)) {
            TreeNode root = tree.getRoot();                    // Retrieves the root,
            node.parent = root;                                 // Parent of this node -> root.
            root.children.add(node);                            // Adds this node as the children of the root.
            tree.nodeList.add(node);                            // Adds this node to the list of nodes.
            assert tree.nodeList.get(0).name.equals("root");
            tree.parentList.put(node.name, 0);                  // Adds "this node -> root" parent mapping.
        }

        // Constructs internal nodes for the rest of the majority bi-partitions and rewires them.

        for (int z = 1; z < tree.clusters.size(); z++) {
            Cluster cluster = tree.clusters.get(z);
            // only take the majority ones....!
            if(cluster.isMajority == true){
	
	            // 1. Retrieves the parent of the first node in this cluster.
	            TreeNode parent = tree.nodeList.get(tree.parentList.get(cluster.get(0).name));
	
	            // 2. Constructs a new internal node.
	            String nodeName = "int" + Integer.toString(numOfNode);
	            TreeNode internalNode = new TreeNode(nodeName);
	            internalNode.addProperty("noOfOccurrences", cluster.noOfOccurrences);	            
	            internalNode.edgeLength = cluster.edgeLength;
	            internalNode.parent = parent;
	
	            // 3. Insert the new node into the node list.
	            tree.nodeList.add(internalNode);
	            assert tree.nodeList.get(tree.nodeList.size() - 1).name.equals(internalNode.name);
	            tree.parentList.put(nodeName, tree.nodeList.size() - 1);
	            
	            // update the clusters node references...
	            // (use a method of storing edges that referenced the positions in the Tree's nodelists for reference later.)
	            tree.clusters.get(z).nodeRefA = tree.parentList.get(cluster.get(0).name);
	            tree.clusters.get(z).nodeRefB = (tree.nodeList.size()-1);        
	            
	            for (TreeNode node : cluster) {
	                // 4. Makes this node the child of the new internal node.
	                node.parent = internalNode;
	                assert node.parent.name.equals(tree.nodeList.get(tree.nodeList.size() - 1).name);
	                tree.parentList.put(node.name, tree.nodeList.size() - 1);
	                internalNode.children.add(node);
	
	                // 5. Delete the moved node(s) from the parent's children.
	                // TODO: optimize? probably not.
	                for (int i = 0; i < parent.children.size(); i++) {
	                    if (parent.children.get(i).name.equals(node.name)) {
	                        parent.children.remove(i);
	                        break;
	                    }
	                }
	            }
	            // Wires up the internal node.
	            parent.children.add(internalNode);
	            numOfNode++;
            }
        }
        return tree;
    }
    
    /**
     * Constructs the network beginning with the consensus tree and then calls the network drawing function to calculate positions
     * 
     * @param tree Consensus tree already calculated
     */
    public CNetwork constructNetwork(CTree tree){
    	// Create a new network...
    	CNetwork network = new CNetwork(noOfTaxa);
        // Let's copy over all of the splits from the tree into the network format.
    	// Until noted, the node lists in tree and network MUST match and the first lot of nodes (after the root) in the tree are assumed to be the taxa.
    	for(int i=0;i<tree.nodeList.size();i++){
    		TreeNode treeNode = tree.nodeList.get(i);
    		// create a new node
    		CNetworkNode networkNode = new CNetworkNode();
    		// copy over the required data and add it...
    		networkNode.Taxaname = treeNode.name;
    		network.nodes.add(networkNode);  		
    	}
    	// Now add in the splits... (ignore 1st cluster as is for star tree)
    	for(int i=1;i<tree.clusters.size();i++){
    		Cluster cluster = tree.clusters.get(i);
    		if (cluster.isMajority == true){
	    		CNetworkSplit networkSplit = new CNetworkSplit();
	    		CNetworkEdge networkEdge = new CNetworkEdge();
	    		networkEdge.split = networkSplit;
	    		networkEdge.networkNodeA = network.nodes.get(cluster.nodeRefA);
	    		networkEdge.networkNodeB = network.nodes.get(cluster.nodeRefB);
	    		networkSplit.edges.add(networkEdge);
	    		// Maybe we shouldn't be copying this over... but it at least clarifies things in code.
	    		networkSplit.edgelength = cluster.edgeLength;
	    		networkSplit.noOfOccurences = cluster.noOfOccurrences;
	    		networkSplit.split = cluster.aboveSplit;
	    		network.nodes.get(cluster.nodeRefA).joins.add(networkEdge);
	    		network.nodes.get(cluster.nodeRefB).joins.add(networkEdge);
	    		network.splits.add(networkSplit);
	    		cluster.added = true;
    		}
    	}
    	// Now add in the taxa length splits.... NOTE i=0 is root?!  Tried and tested anyway...
    	for(int i=1;i<=noOfTaxa;i++){
    		CNetworkSplit networkSplit = new CNetworkSplit();
    		CNetworkEdge networkEdge = new CNetworkEdge();
    		networkEdge.split = networkSplit;
    		int taxonParent = tree.parentList.get(tree.nodeList.get(i).name);
    		networkEdge.networkNodeA = network.nodes.get(taxonParent);
    		networkEdge.networkNodeB = network.nodes.get(i);
    		networkSplit.edges.add(networkEdge);
    		networkSplit.edgelength = tree.nodeList.get(i).edgeLength;
    		networkSplit.noOfOccurences = noOfTrees;
    		BitSet trivialSplit = new BitSet(noOfTaxa);
    		trivialSplit.set(i-1);
    		networkSplit.split = trivialSplit;
    		network.nodes.get(taxonParent).joins.add(networkEdge);
    		network.nodes.get(i).joins.add(networkEdge);
    		network.splits.add(networkSplit);
    	}
    	// NOTE: From now on the indices on Network & Tree may no longer match!!
    	// Remove the root if it is simply disecting an edge (only 2 edges on either side).
		// NOTE: root is assumed to be the first position in the node list in CNetwork and CTree
    	if(network.nodes.get(0).joins.size() == 2){
    		// find the nodes on either side of the root...
    		CNetworkNode[] redundant = new CNetworkNode[2];
    		for(int h=0;h<2;h++){
				if(network.nodes.get(0).joins.get(h).networkNodeA == network.nodes.get(0)){
					redundant[h] = network.nodes.get(0).joins.get(h).networkNodeB;
				}
				else{
					redundant[h] = network.nodes.get(0).joins.get(h).networkNodeA;
				}
    		}
    		// correct the edge in the one list to go to the zero node.
    		if(network.nodes.get(0).joins.get(1).networkNodeA == network.nodes.get(0)){
    			network.nodes.get(0).joins.get(1).networkNodeA = redundant[0];
    		}
    		else{
    			network.nodes.get(0).joins.get(1).networkNodeB = redundant[0];
    		}
    		// find the edge in the list in zero and remove it..
    		for(int h=0;h<redundant[0].joins.size();h++){
    			if(redundant[0].joins.get(h).networkNodeA == network.nodes.get(0) || redundant[0].joins.get(h).networkNodeB == network.nodes.get(0)){
    				redundant[0].joins.remove(h);
    				h--;
    			}
    		}
    		// remove the edge from the split list
    		for(int h=0;h<network.nodes.get(0).joins.get(0).split.edges.size();h++){
    			if(network.nodes.get(0).joins.get(0).split.edges.get(h) == network.nodes.get(0).joins.get(0)){
    				network.nodes.get(0).joins.get(0).split.edges.remove(h);
    				h--;
    			}
    		}
    		// add the other edge in the zero list
    		redundant[0].joins.add(network.nodes.get(0).joins.get(1));
    		// remove the node at last, and any reference to the removed edge also goes!
    		network.nodes.remove(0);
    	}

    	// TODO: Optional: sort the list in terms of frequency of occurrence to select only most frequent if less than 33.34% occurences....!
        //IGNORE the first one as it is the star network...
    	//Now let's go through the remaining splits and try and add them into the network.
    	for(int i=0;i<tree.clusters.size();i++){
    		Cluster cluster = tree.clusters.get(i);
    		// for clusters that haven't been added yet, let's add them...
    		if (cluster.added == false){
    			//define new split and copy over info.
	    		CNetworkSplit addSplit = new CNetworkSplit(); 
	    		addSplit.split = cluster.aboveSplit;
	    		addSplit.edgelength = cluster.edgeLength;
	    		// now findout if is incompatible with any existing splits in the network
	    		ArrayList<CNetworkSplit> incompatibleSplits = new ArrayList<CNetworkSplit>();
	    		for(CNetworkSplit compSplit : network.splits){
		    		if(compSplit.isCompatible(addSplit,noOfTaxa)==false){
		    			incompatibleSplits.add(compSplit);
		    		} 			
	    		}
	    		// Now we need to find the nodes we need to copy, so create a list of them
	    		ArrayList<CNetworkNode> nodesToCopy = new ArrayList<CNetworkNode>();
    			// go through the incompatible splits one by one and add them to chains.
    			ArrayList<CNetworkPath> chains = new ArrayList<CNetworkPath>();
    			// note all nodes that will be stored in paths for more efficient path claering up
    			ArrayList<CNetworkNode> nodesToClearUp = new ArrayList<CNetworkNode>();
    			// Loop through incompatible splits
    			for(CNetworkSplit currentSplit : incompatibleSplits){
    				// Loop through all the edges in this split
    				for(CNetworkEdge currentEdge : currentSplit.edges){
    					// create a new path and add it to the list...
						CNetworkPath newPath = new CNetworkPath();
						newPath.pathOfEdges.add(currentEdge);
						newPath.pathOfNodes.add(currentEdge.networkNodeA);
						newPath.pathOfNodes.add(currentEdge.networkNodeB);
						currentEdge.networkNodeA.AddToPaths(newPath,nodesToClearUp);
						currentEdge.networkNodeB.AddToPaths(newPath,nodesToClearUp);
						chains.add(newPath);
						// Now extend any existing paths...
    					// create a list of paths to extend found paths with...
    	    			ArrayList<CNetworkPath> extendPathsWith = new ArrayList<CNetworkPath>();
    	    			extendPathsWith.add(newPath);
    					// Hack to go from NodeA & NodeB to a for loop...
    					CNetworkNode[] edgeNodes = new CNetworkNode[]{currentEdge.networkNodeA, currentEdge.networkNodeB};
    					// For the first node, add the standard
    					 for(int l=0; l<=1;l++){
							// get ready to make this part of the nodes with paths to clear up... (should be already done?)
							if(edgeNodes[l].paths.isEmpty()==true){
								nodesToClearUp.add(edgeNodes[l]);
							}
    						// Look at each path that the node is involved in...
							int pathCount = edgeNodes[l].paths.size();
	    					for(int g= 0; g <pathCount; g++){
	    						 CNetworkPath currentPath = edgeNodes[l].paths.get(g);
	    						 // for the case that we may want to add the edge to the start/end of a path on the node...
								 if(currentPath.endNode()==edgeNodes[l] || currentPath.startNode()==edgeNodes[l]){
									 //check we haven't already added this split yet...
									 if (currentPath.checkForSplit(currentEdge)==false){
										 //create a new combined path with all paths to extend with...
										 //Don't forget to make this add the new path to all the nodes as well as combining over the correct node & removing one instance of it...
										for(CNetworkPath pathToAdd : extendPathsWith){
											chains.add(currentPath.combineWith(pathToAdd,edgeNodes[l]));
										}
									 }
								 }
	    					 }
    					 }
    				}
    			}
    			for(CNetworkNode nodeToClearUp : nodesToClearUp){
    				nodeToClearUp.paths.clear();
    			}
    			ArrayList<CNetworkEdge> incompatibleEdges = new ArrayList<CNetworkEdge>();
    			//Case that there are no incompatible splits... just add in the split!
    			// Note that SOME splits appear to contain the whole tree (perhaps due to a hack earlier), so lets just ignore these....  This could be described as a hack induced hack...
    			if(incompatibleSplits.isEmpty() == true && (addSplit.split.cardinality() != noOfTaxa && addSplit.split.cardinality() !=0)){
    				// we have no idea which taxon is actually on the side we want, but we will use some cool functions of BitSet to help reduce the time it takes...
    				CNetworkNode beginNode;
    				// Start searching from the adding side of the split containing zeros?
    				boolean zeroSide;    				
    				if (addSplit.split.cardinality() <= (noOfTaxa)/2){
    					// if there are more set to 0 than 1 go from a taxon on the 1 side:
    					// Search through the nodes for this using predefined function...
    					beginNode = TaxonRefToNode(addSplit.split.nextSetBit(0), network);
    					zeroSide = false;
    				}
    				else{
    					// if there are more set to 1 than 0, then go from a taxon on the 0 side:
    					beginNode = TaxonRefToNode(addSplit.split.nextClearBit(0), network);
    					zeroSide = true;
    				}
    				// recursively go along the edges until a node is found that has an edge that satisfies.
    				ArrayList<CNetworkNode> consideredNodes = new ArrayList<CNetworkNode>();
    				CNetworkNode CopyNode2 = beginNode.findFirstNonSubset(addSplit,consideredNodes, zeroSide,noOfTaxa);
    				nodesToCopy.add(CopyNode2);
    			}
    			else{
    				//Copy over the list of nodes from the longest path...ERROR if we have a too big path...
    				if(chains.size()>0){
    					// Store the max paths...
    					ArrayList<CNetworkPath> currentMaxPaths = new ArrayList<CNetworkPath>();
    					// initialise with first one, which only serves as to store length...:
    							currentMaxPaths.add(chains.get(0));
	    				for(CNetworkPath chain : chains){
	    					if(chain.pathOfNodes.size()>currentMaxPaths.get(0).pathOfNodes.size()){
	    						currentMaxPaths.clear();
	    						currentMaxPaths.add(chain);
	    					}
	    					// Yes, still add it if we already added it!
	    					if(chain.pathOfNodes.size()==currentMaxPaths.get(0).pathOfNodes.size()){
	    						currentMaxPaths.add(chain);
	    					}
	    				}
	    				// now remove any max path that requires all edges to be copied or none at all...
	    				for(int m=1;m<currentMaxPaths.size();m++){
	    					if(currentMaxPaths.get(m).isSplitSide(addSplit, noOfTaxa)==false){
	    						currentMaxPaths.remove(m);
		    					m--;
	    					}
	    				}
	    				if(currentMaxPaths.size()>2){
	        				//ERROR!!! Now in >2d....!
	    					System.out.println("Potential Error: Network now in > 2 dimensions....!");
	        			}
	    				// Now create the copyNode list:
	    				// Store those that have been copied.... err
	    				// Go through all the paths
	    				// DOES START AT 1 as first just stored length!
	    				for(int m=1;m<currentMaxPaths.size();m++){
	    					incompatibleEdges.addAll(currentMaxPaths.get(m).pathOfEdges);
	    					CNetworkPath currentMaxPath = currentMaxPaths.get(m);
	    					// Go through all the nodes
	    					for(CNetworkNode currentNodeCopy : currentMaxPath.pathOfNodes){
	    						if(currentNodeCopy.consideration == false){
	    							currentNodeCopy.consideration = true;
	    							nodesToCopy.add(currentNodeCopy);
	    						}
	    					}
	    				}
	    				// remove the considered flag from the node:
	    				for(CNetworkNode currentNodeClean : nodesToCopy){
	    					currentNodeClean.consideration = false;
	    				}
    				}

    			}
    			// actually call the function that copies the nodes
    			CopyNode(nodesToCopy,addSplit,network,incompatibleEdges);
    			// and add the new split at last to the network.
    			network.splits.add(addSplit);
    		}
    	}
        // fill in the network positions (hopefully)...
        network.FindPositions(noOfTaxa);
        //printDetails(tree,network);
    	return network;
    }

    
    /**
     * Copy the list of nodes that need copying when adding a new split into a network.
     *  
     * @param nodesToCopy ArrayList of nodes in the network that need to be copied
     * @param network Network we are working in
     * @param incompatibleEdges Edges that are in the path that joins up the nodes that are being copied
     * @param splitToAdd Split that we are currently adding in by copying these nodes
     *
     */
    public void CopyNode(ArrayList<CNetworkNode> nodesToCopy, CNetworkSplit splitToAdd, CNetwork network, ArrayList<CNetworkEdge> incompatibleEdges){
    	// Copy Node:
    	for(CNetworkNode nodeToCopy: nodesToCopy){
    		CNetworkNode newNode = new CNetworkNode();
    		// create a new node with a name that indicates it has been added (we can output node names by modifying later code)
    		newNode.Taxaname = "added0"+Integer.toString(network.nodes.size());
    		// below should NEVER happen... but will leave in just in case as was at one time an issue:
	    	if (nodeToCopy == null || newNode == null){
	    		System.out.println("Null pointer in terms of node that we want to copy to...");
	    	}
	    	// Reference the new node on the old one
    		nodeToCopy.CopiedToNode = newNode;
    		// Add new node to network
    		network.nodes.add(newNode);
    		// Add an edge between old and new node
    		CNetworkEdge edgeToAdd = new CNetworkEdge();
    		edgeToAdd.networkNodeA = newNode;
    		edgeToAdd.networkNodeB = nodeToCopy;
    		edgeToAdd.split = splitToAdd;
    		splitToAdd.edges.add(edgeToAdd);
    	}
    	// Go through and move links over depending on which side the splits are a subset of.
    	for(CNetworkNode nodeToCopy: nodesToCopy){
    		// look at each join on each old node
			for(int i=0;i<nodeToCopy.joins.size();i++){
    			CNetworkEdge edgeToConsider = nodeToCopy.joins.get(i);
    			// if it joins another copied node then we will duplicate the edge later...
    			if(edgeToConsider.networkNodeA.CopiedToNode == null || edgeToConsider.networkNodeB.CopiedToNode == null){
            		//if on "0" side of new split then keep it, otherwise move it to the new node.
    				if(splitToAdd.OnZeroSide(edgeToConsider.split,noOfTaxa)==false){
    					// Replace the appropriate edge reference...
    					if(edgeToConsider.networkNodeA == nodeToCopy){
    						edgeToConsider.networkNodeA = nodeToCopy.CopiedToNode;
    						//System.out.println("moved link to "+edgeToConsider.networkNodeB.Taxaname+" from "+nodeToCopy.Taxaname+" to "+nodeToCopy.CopiedToNode.Taxaname);
    					}
    					else if(edgeToConsider.networkNodeB == nodeToCopy){
    						edgeToConsider.networkNodeB = nodeToCopy.CopiedToNode;
    						//System.out.println("moved link to "+edgeToConsider.networkNodeA.Taxaname+" from "+nodeToCopy.Taxaname+" to "+nodeToCopy.CopiedToNode.Taxaname);
    					}
    					// Add it to the new node list
						nodeToCopy.CopiedToNode.joins.add(edgeToConsider);
    					// delete it from the list of the old nodes:
						nodeToCopy.joins.remove(i);
						// go back one in the list!
						i--;
    				}
    			}
    		}
    	}
    	// Add in the edges between copied nodes to the old nodes to their node join lists now we have done all the moving about!
    	for(CNetworkEdge edgeToLink : splitToAdd.edges){
    		edgeToLink.networkNodeA.joins.add(edgeToLink);
    		edgeToLink.networkNodeB.joins.add(edgeToLink);
    	}
    	// Add in the edges to the copied nodes that are in the list of incompatible splits.
		for(CNetworkEdge incompatibleEdge : incompatibleEdges){
			// Check to see if we have an edge that is incompatible (this SHOULD always be the case... )
			if(incompatibleEdge.networkNodeA.CopiedToNode !=null && incompatibleEdge.networkNodeB.CopiedToNode != null){
				// Create new edge
				CNetworkEdge newEdge = new CNetworkEdge();
				// Update node info.
				newEdge.networkNodeA = incompatibleEdge.networkNodeA.CopiedToNode;
				newEdge.networkNodeB = incompatibleEdge.networkNodeB.CopiedToNode;
				newEdge.networkNodeA.joins.add(newEdge);
				newEdge.networkNodeB.joins.add(newEdge);
				// Update split info.
				newEdge.split = incompatibleEdge.split;
				incompatibleEdge.split.edges.add(newEdge);
			}
		}
    	// Now remove the CopiedToNode link now we are done...
    	for(CNetworkNode nodeToCopy: nodesToCopy){
    		nodeToCopy.CopiedToNode = null;
    	}
    	
    }
    /**
     * Find the taxon node in the network.
     *  
     * @param network The network to find the taxon in
     * @param taxonRef Integer representation of taxon in the TaxaMap
     *
     */
    public CNetworkNode TaxonRefToNode(int taxonRef, CNetwork network){
    	//function to return the node that is of the taxon
    	for(CNetworkNode currentNode : network.nodes){
    		if(currentNode.Taxaname == taxa.getName(taxonRef)){
    			return currentNode;
    		}
    	}
    	// Notify that method failed... return null...
    	System.out.println("Error - taxon reference not found.  Maybe duplicate names?");
    	return null;
    }

    public static String[] readTreesFromFile(String fileName, int n) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fileName));
        String[] trees = new String[n];
        for (int i = 0; i < n; i++) {
            trees[i] = reader.readLine();
        }
        return trees;
    }

    private static void printUsage() {
        System.out.println("Usage: java CTMain <file containing trees> <no. of trees>");
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 2) {
            printUsage();
            System.exit(-1);
        }

        //CTMain main = new CTMain();
        //String[] trees = readTreesFromFile(args[0], Integer.parseInt(args[1]));
        //System.out.println(main.start(trees));
    }

    // Getters and setters

    public double getResPercentage() {
        return resPercentage;
    }

    public void setResPercentage(double resPercentage) {
        this.resPercentage = resPercentage;
    }

    public long getSeed() {
        return seed;
    }

    public void setSeed(long seed) {
        this.seed = seed;
    }

    //THE following 3 functions are used for testing only....
    /**
     * FOR TESTING ONLY
     * Function called to initialise the Network Tester hash table & rest of CTMain
     * 
     * @param noOfTestTaxa Number of taxa
     * @param noOfTestSamples number of trees input
     *
     */
    public void InitialiseNetworkTester(int noOfTestTaxa,int noOfTestSamples){
        // TaxaMap initialisation
        noOfTrees = 0;
        noOfTaxa = noOfTestTaxa;
        noOfSamples = noOfTestSamples;
        // Hash initialisation
        hashUtils = new HashUtils();
        hashUtils.initialize(noOfTaxa, noOfSamples, C, seed);
        hashTable = new HashTable(hashUtils.m1);
        leafEdgeLengths = new double[noOfTaxa];
        taxa = new TaxaMap(noOfTaxa);
        for (int i = 0; i < noOfTaxa; i++) {
            taxa.put(""+i, i);
            leafEdgeLengths[i] = 1;
        }


        // Adds a single star partition, once and for all.
        BitSet star = new BitSet(noOfTaxa);
        star.flip(0, noOfTaxa);
        HashEntry entry = new HashEntry(-1, star, 0.0d);
        entry.count = noOfSamples+1;
        partitions.add(entry);
        updateInterestThreshold(); 
	}
    /**
     * FOR TESTING ONLY
     * Add a tree to the hash table
     * 
     * @param splits Tree to add consisting of an arraylist of BitSet splits
     *
     */
    public void addTestTree(ArrayList<BitSet> splits){
        noOfTrees++;
        updateInterestThreshold();
        // Hashes the partitions of the trees.
        noOfPartitions = 0;
        for(BitSet partition : splits){
        	addTestSplit(partition);
        }
	}
    /**
     * FOR TESTING ONLY
     * Add a single split to the hash table
     * 
     * @param partition BitSet split to add
     *
     */
	public void addTestSplit(BitSet partition){
		long tableKey2 = 0; 
        long bucketKey2 = 0; 
        // always calculate the split directly
		for (int k=0;k<noOfTaxa;k++) { 
	        if(partition.get(k)==true){ 
	                tableKey2 += hashUtils.a1[k]; 
	                bucketKey2 += hashUtils.a2[k]; 
	        } 
		} 
		hashTable.put(partition, 1, (int) (tableKey2 % hashUtils.m1), (int) (bucketKey2 % hashUtils.m2), interestThreshold,partitions); 
	}
	
    /**
     * FOR TESTING ONLY
     * Print details from the network - call after it has been constructed!
     * 
     * @param tree Standard consensus tree
     * @param network The network
     *
     */
	public void printDetails(CTree tree,CNetwork network){
    	// Print out the clusters:
       	System.out.println("ALL CLUSTERS:");
    	for(int i=0;i<tree.clusters.size();i++){
    		System.out.println(i+" "+tree.clusters.get(i).aboveSplit.toString()+" "+tree.clusters.get(i).noOfOccurrences);
    	}
    	// Print out the tree copied into network format:
    	System.out.println("TREE IN NETWORK, splits:");
    	for(int i=0;i<network.splits.size();i++){
    		System.out.println(i+" "+network.splits.get(i).split.toString()+" "+network.splits.get(i).noOfOccurences);
        	for(int j=0;j<network.splits.get(i).edges.size();j++){
        		System.out.println(j+" "+network.splits.get(i).edges.get(j).networkNodeA.Taxaname+" to "+network.splits.get(i).edges.get(j).networkNodeB.Taxaname);
        	}
    	}
    	 System.out.println("The following comes from the network format...");
         System.out.println("NODES:");
         for (int i = 0; i< network.nodes.size();i++){
             System.out.println(""+network.nodes.get(i).hashCode()+" "+network.nodes.get(i).Taxaname);
         }
         
         System.out.println("EDGES:");
         for (int i = 0; i < network.splits.size(); i++){
         	for(int j=0; j<network.splits.get(i).edges.size();j++){
         		System.out.println(network.splits.get(i).split.toString()+": "
         	+network.splits.get(i).edges.get(j).networkNodeA.Taxaname+" TO "+
         	network.splits.get(i).edges.get(j).networkNodeB.Taxaname);
         	}
         }
         System.out.println("The following comes from the original tree...(!)");
         for (int i = 1; i < tree.clusters.size(); i++) {
         	Cluster cluster2 = tree.clusters.get(i);
             //if(cluster2.isMajority == false){
 	        	network.outputString = "1 ";
 	        	for(TreeNode node :cluster2){
 	        		network.outputString += " " + (taxa.get(node.name)+1);
 	        	}
 	    		network.outputString += ",";
 	            System.out.println(network.outputString);
             //}
         }
       //printing to a file
         try {
 			network.PrintOut(1,noOfTaxa,taxa);
 		} catch (IOException e) {
 			// TODO Auto-generated catch block
 			e.printStackTrace();
 		}
         //network in nex format
         try {
         	network.PrintOut(2,noOfTaxa,taxa);
 		} catch (IOException e) {
 			// TODO Auto-generated catch block
 			e.printStackTrace();
 		}
    	return;
	}

}
