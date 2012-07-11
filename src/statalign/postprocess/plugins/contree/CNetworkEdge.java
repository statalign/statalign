package statalign.postprocess.plugins.contree;

import java.util.ArrayList;
/**
 * Edges in the consensus network
 * 
 * @author wood
 *
 */
public class CNetworkEdge {
	//store the nodes that it starts and ends at.
	public CNetworkNode networkNodeA;
	public CNetworkNode networkNodeB;
	// the corresponding Split this edge is part of
	public CNetworkSplit split;
	// list of neighbours when calculating positions (i.e. edges it will form boxes with)
	ArrayList<CNetworkEdge> neighbours;
	// booelan for consideration in drawing algorithm at a node when assigning to neighbour lists...
	public boolean considered;
	// weighting of edge when drawing determining how much angle it will receive
	public double weighting;
	// end positions
	public double xPosA,yPosA,xPosB,yPosB;
	// boolean true if drawn (has x, y positions assigned), otherwise false
	public boolean drawn;
	/**
	 * Set up a new edge
	 *
	 */
	public CNetworkEdge(){
		drawn = false;
		considered = false;
		xPosA = -1.0;
		xPosB = -1.0;
		yPosA = -1.0;
		yPosB = -1.0;		
	}
}