package statalign.postprocess.plugins.contree;

import java.util.ArrayList;
/**
 * Paths of incompatible splits in the consensus network
 * 
 * @author wood
 *
 */
public class CNetworkPath {
	public ArrayList<CNetworkEdge> pathOfEdges;
	public ArrayList<CNetworkNode> pathOfNodes;
	public CNetworkPath(){
		pathOfEdges = new ArrayList<CNetworkEdge>();
		pathOfNodes = new ArrayList<CNetworkNode>();
	}
	/**
	 * Returns the node at the start of the path
	 */
	public CNetworkNode startNode(){
		return pathOfNodes.get(0);
	}
	/**
	 * Returns the node at the end of the path
	 */
	public CNetworkNode endNode(){
		return pathOfNodes.get(pathOfNodes.size()-1);
	}
	/**
	 * Returns true if this path already contains an edge from the argument split, otherwise false.
	 * 
	 * @param edgeToLookFor Edge from split we want to see is already present in path or not
	 */
	public boolean checkForSplit(CNetworkEdge edgeToLookFor){
		for(CNetworkEdge currentEdge : pathOfEdges){
			if(currentEdge.split == edgeToLookFor.split){
				return true;
			}
		}
		return false;
	}
	/**
	 * Combines this path with the argument path across the argument node
	 * 
	 * @param toAdd Path to add to this one
	 * @param nodeToLinkOver Node over which to link (ensures non-duplication of this node)
	 */
	public CNetworkPath combineWith(CNetworkPath toAdd , CNetworkNode nodeToLinkOver){
		// Path to Return (does what it says...)
		CNetworkPath pathToReturn = new CNetworkPath();
		//we'll always add to this one
		if(pathOfNodes.get(0)== nodeToLinkOver){
			CNetworkPath firstPath = reordered();
			pathToReturn.pathOfEdges.addAll(firstPath.pathOfEdges);
			pathToReturn.pathOfNodes.addAll(firstPath.pathOfNodes);
		}
		else{
			pathToReturn.pathOfEdges.addAll(pathOfEdges);
			pathToReturn.pathOfNodes.addAll(pathOfNodes);
		}
		// now we have the path with the node to link over at the end.
		// and now add the second one given as an argument
		// note we are careful not to modify the already existing paths as we want more paths.
		if(toAdd.pathOfNodes.get(0) == nodeToLinkOver){
			for(int j=1;j<toAdd.pathOfNodes.size();j++){
				pathToReturn.pathOfNodes.add(toAdd.pathOfNodes.get(j));
			}
			pathToReturn.pathOfEdges.addAll(toAdd.pathOfEdges);
		}
		else{
			CNetworkPath secondPath = toAdd.reordered();
			secondPath.pathOfNodes.remove(0);
			pathToReturn.pathOfEdges.addAll(secondPath.pathOfEdges);
			pathToReturn.pathOfNodes.addAll(secondPath.pathOfNodes);

		}
		// add to all the nodes in the path this path!
		for(CNetworkNode nodeToAddPathTo : pathOfNodes){
			nodeToAddPathTo.paths.add(pathToReturn);
		}
		return pathToReturn; 
	}
	/**
	 * Reorders this path
	 */
	public CNetworkPath reordered(){
		CNetworkPath reorderedPath = new CNetworkPath();
		for(int i=pathOfEdges.size() - 1;i>=0;i--){
			reorderedPath.pathOfEdges.add(pathOfEdges.get(i));
		}
		for(int i=pathOfNodes.size() - 1;i>=0;i--){
			reorderedPath.pathOfNodes.add(pathOfNodes.get(i));
		}
		return reorderedPath;
	}
	/**
	 * This function sees if adjoining splits at the nodes at the path ends create subsets on both sides of the input split...
	 * Is to see is path actually links two nodes from where the split occurs.
	 * Returns true if is a desired split path otherwise false.
	 * 
	 * @param inputSplit Split that is being added and must create splits in some of the edges joining the end nodes
	 * @param noOfTaxa Standard request for number of taxa in the network
	 */
	// 
	public boolean isSplitSide(CNetworkSplit inputSplit, int noOfTaxa){
		// booleans to store if we have edges with splits as subsets of zeroside and one side of the input split 
		boolean zeroSideSplit = false;
		boolean oneSideSplit = false;
		// look at both ends of the path (code hack...)
		CNetworkNode[] nodesToCheck = new CNetworkNode[]{pathOfNodes.get(0),pathOfNodes.get(pathOfNodes.size()-1)};
		for(int l=0;l<=1;l++){
			// look at all joins...
			for(CNetworkEdge currentEdge : nodesToCheck[l].joins){
				// ignore incompatible splits in the path... 
				if(checkForSplit(currentEdge) == false){
					// and now decide which side it is on
					if(inputSplit.isSubset(currentEdge.split, true, noOfTaxa)){
						zeroSideSplit = true;
					}
					else if(inputSplit.isSubset(currentEdge.split, false, noOfTaxa)){
						oneSideSplit = true;
					}
				}
			}	
		}
		// Return true if is split path
		if(zeroSideSplit==true && oneSideSplit==true){
			return true;
		}
		// otherwise false!
		else{
			return false;
		}
	}
}