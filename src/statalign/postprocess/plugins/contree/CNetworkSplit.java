package statalign.postprocess.plugins.contree;

import java.util.ArrayList;
import java.util.BitSet;
/**
 * Splits in the consensus network
 * Includes several functions to check degrees of compatibility with another split
 * 
 * @author wood
 *
 */
public class CNetworkSplit {
	// Edges that represent this split
	public ArrayList<CNetworkEdge> edges;
	// The split in binary format
	public BitSet split;
	// The length of this split
	public double edgelength;
	// The weight based on the proportion of times that this edge occurs in the tree samples.
	public int noOfOccurences;
	// Direction of split in angle, is negative if unset.
	public double direction = -1.0;
	/**
	 * Constructs a new split.
	 */
	public CNetworkSplit(){
		edges = new ArrayList<CNetworkEdge>();
	}
	/**
	 * Compares an input split to see if it is compatible with this one.  i.e. is one of the sides of one split entirely within one side of the other.
	 * 
	 * @param comp The split to compare this one with
	 * @param noOfTaxa The number of taxa as the count function on size of a BitSet is strange
	 */
	public boolean isCompatible(CNetworkSplit comp, int noOfTaxa){
		BitSet split1 = this.split;
		BitSet split2 = comp.split;
		// See if any side can be considered a subset of another.  If so, then is compatible.
		// Do this by seeing if a 1 in split1 is always either a 1 or zero in split2 OR
		// a 0 in split1 is always either a 1 or zero in split2
		boolean[] combinations = new boolean[]{false,false,false,false};
		for(int i=0; i<noOfTaxa;i++){
			combinations[((split1.get(i)==true)?1:0)*2+((split2.get(i)==true)?1:0)] = true;
		}
		if (combinations[0] == false || combinations[1] == false || combinations[2] == false || combinations[3] == false){
			return true;
		}
		else{
			return false;
		}
	}
	/**
	 * Compares an input split to see if the input split occurs on the zero side of this split.
	 * 
	 * @param comp The split to try and find in the zero side of this split
	 * @param noOfTaxa The number of taxa as the count function on size of a BitSet is strange
	 */
	public boolean OnZeroSide(CNetworkSplit comp, int noOfTaxa){
		BitSet split1 = this.split;
		BitSet split2 = comp.split;
		// We only want to return TRUE if comparison split is a subset of the zero side of the split calling the function (This split)#
		// We can test for this by seeing if all the "1"s in this split match with a single value (of either "0" or "1" in the comparison
		boolean[] found = new boolean[]{false,false};
		for(int i=0; i<noOfTaxa;i++){
			if(split1.get(i)==true){
				found[(split2.get(i)==true)?1:0] = true;
			}
		}
		if(found[0] == true && found[1] == true) return false;
		else return true;
	}
	/**
	 * Compares an input split to see if the argument split is a subset of the specified side of this split?
	 * NOTE: assumes that splits are COMPATIBLE to begin with i.e. all combinations CANNOT be true!!
	 * 
	 * @param comp The split we want to find as a subset
	 * @param zeroSide If we want the split on the zeroside of this split then set as true, if on the other side then as false...
	 * @param noOfTaxa The number of taxa as the count function on size of a BitSet is strange
	 */
	public boolean isSubset(CNetworkSplit comp, boolean zeroSide, int noOfTaxa){
		BitSet split1 = this.split;
		BitSet split2 = comp.split;
		// shamelessly copy code from the isCompatible function...
		boolean[] combinations = new boolean[]{false,false,false,false};
		for(int i=0; i<noOfTaxa;i++){
			combinations[((split1.get(i)==true)?1:0)*2+((split2.get(i)==true)?1:0)] = true;
		}
		// for the case that we have all the splits on one side
		if(zeroSide == true && combinations[0] == true && combinations[1] == true){
			return true;
		}
		else if (zeroSide == false && combinations[2] == true && combinations[3] == true){
			return true;
		}
		else return false;
	}
	/**
	 * Takes an input edge and finds the side of this split it occurs on and then returns the number of taxa on that side of this split. 
	 * NOTE: assumes that splits are COMPATIBLE to begin with i.e. all combinations CANNOT be true!
	 * 
	 * @param currentEdge An edge from this node to compare
	 * @param noOfTaxa The number of taxa as the count function on size of a BitSet is strange
	 */
	// Assumed that splits input are compatible...
	public double getSubsetJoins(CNetworkEdge currentEdge,int noOfTaxa){
		CNetworkSplit splitToCompare = currentEdge.split;
		// if is a subset on ZERO side, then return number of zeros...
		if(isSubset(splitToCompare,true,noOfTaxa)==true) {
			return (double)(noOfTaxa-split.cardinality());
		}
		// if on 1 side, return the number of 1s....
		else if (isSubset(splitToCompare,false,noOfTaxa)==true){
			return (double)split.cardinality();
		}
		// shouldn't be returned... but will be if, for example, a node in the middle of a line is output.
		else{
			System.out.println("Potential problem in drawing network: No subset found for splits "+splitToCompare.split.toString()+" and "+split.toString());
			return currentEdge.weighting;
		}
	}
}