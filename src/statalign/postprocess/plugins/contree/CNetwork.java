

package statalign.postprocess.plugins.contree;



import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * The Consensus network
 * 
 * @author wood
 *
 */
public class CNetwork {
	public String outputString;
	public ArrayList<CNetworkSplit> splits;
	public ArrayList<CNetworkNode> nodes;
	public int noOfTaxa;
	/**
	 * Constructs a new Consensus Network.
	 * 
	 * @param noOfTaxa The number of taxa in the network.
	 */
	public CNetwork(int noOfTaxa){
		splits = new ArrayList<CNetworkSplit>();
		nodes = new ArrayList<CNetworkNode>();
		this.noOfTaxa = noOfTaxa;
	}
	/**
	 * Prints out the network in Nexus format.  Can be read by SplitsTree 4
	 * 
	 * @param function 1 - print as splits, 2 - print drawn network as calculated
	 * @param noOfTaxa The number of taxa in the network.
	 * @param taxa the map of taxa indices to names
	 */
    public void PrintOut(int function,int noOfTaxa,TaxaMap taxa) throws IOException{
    	String fileName = "";
    	 if(function==1){
    		 fileName = "splits";
    	 }
    	 if(function==2){
    		 fileName = "network";
    	 }
        PrintWriter out = new PrintWriter(new FileWriter(fileName+".nex")); 
        out.println("#nexus"); 
        out.println("BEGIN Taxa;");
        out.println("DIMENSIONS ntax="+noOfTaxa+";");
        out.println("TAXLABELS");
        for(int x=0;x<taxa.size();x++){
        	out.println("["+(x+1)+"] '"+taxa.getName(x)+"'");
        }
        out.println(";");
        out.println("END; [Taxa]");
        if(function==1){
	        out.println("BEGIN Splits;");
	        int countSplits = 0;
	        for(int y=0; y<splits.size();y++){
	        	if(splits.get(y).split.cardinality()!=0 && splits.get(y).split.cardinality()!=noOfTaxa){
	        		countSplits++;
	        	}
	        }
	        out.println("DIMENSIONS ntax="+noOfTaxa+" nsplits="+countSplits+";");
	        
	        out.println("FORMAT labels=no weights=yes confidences=no intervals=no;");
	        out.println("PROPERTIES fit=-1.0 compatible;");
	        out.print("CYCLE");
	        for(int x=0;x<taxa.size();x++){
	        	out.print(" "+(x+1));
	        }
	        out.println(";");
	        out.println("MATRIX");
	
	        for(int y=0; y<splits.size();y++){
	        	if(splits.get(y).split.cardinality()!=0 && splits.get(y).split.cardinality()!=noOfTaxa){
		        	out.print("["+(y+1)+", size=1]  	 1.0 	  ");
		        	for(int z=0;z<noOfTaxa;z++){
		        		out.print((splits.get(y).split.get(z))?(z+1)+" ":"");
		        	}
		            out.println(",");
	        	}
	        }
	        out.println(";");
	        out.println("END; [Splits]");
        }
        if(function==2){
        	out.println("BEGIN Network;");
        	int noOfEdges = 0;
        	for(int w=0;w<splits.size();w++){
        		for(int u=0;u<splits.get(w).edges.size();u++){
        			noOfEdges++;
        		}
        	}
        	out.println("DIMENSIONS ntax="+noOfTaxa+" nvertices="+nodes.size()+" nedges="+noOfEdges+";");
        	out.println("DRAW to_scale;");
        	out.println("VERTICES");
        	for(int w=0;w<nodes.size();w++){
            	out.println((w+1)+" "+nodes.get(w).xPos+" "+nodes.get(w).yPos+",");	
            	nodes.get(w).outputNumber = (w+1); 	
        	}
            out.println(";");
        	out.println("VLABELS");
        	for(int w=0;w<nodes.size();w++){
            	out.println((w+1)+" '"+nodes.get(w).Taxaname+"',");	

        	}
        	out.println(";");
        	out.println("EDGES");
        	int edgePrinted = 0;
        	for(int w=0;w<splits.size();w++){
        		for(int u=0;u<splits.get(w).edges.size();u++){
        			edgePrinted++;
                	out.println(edgePrinted+" "+splits.get(w).edges.get(u).networkNodeA.outputNumber+" "+splits.get(w).edges.get(u).networkNodeB.outputNumber+",");
        		}
        	}
        	out.println(";");
        	out.println("END; [Network]");

        }
        out.close();
    	
    }
	/**
	 * Calculates positions for network nodes and edges linking them to output to the network visualiser.
	 * 
	 * @param noOfTaxa The number of taxa in the network
	 */
    public void FindPositions(int noOfTaxa){
    	//catch empty network...
    	if (nodes.size()>0){
    		//pick a node to start with... (we'll just start from the end as this is most likely to be not a taxon (unless it is a star graph hmmmm).
	    	CNetworkNode firstNode = nodes.get(nodes.size()-1);
	    	// it has all angles to go at...
	    	firstNode.angleWidth = 2*Math.PI;
	    	firstNode.angleDirection = 0.0;
	    	firstNode.drawn = true;
	    	// begin at origin...
	    	firstNode.xPos = 0.0;
	    	firstNode.yPos = 0.0;
	    	// the next edges are stored in an ArrayList so that the drawing spreads out from our initial node!
	    	ArrayList<CNetworkNode> nextNodesToDraw = new ArrayList<CNetworkNode>();
	    	nextNodesToDraw.add(firstNode);
	    	// Loop through nodes in rings that spread outwards from the central node, with each ring containing nodes that are the same number of nodes away from the initial node.
	    	do{
	    		// Move next lot over to current lot to draw
	    		ArrayList<CNetworkNode> nodesToDraw = new ArrayList<CNetworkNode>();
	    		nodesToDraw.addAll(nextNodesToDraw);
	    		// and clear the array that stores the next lot
	    		nextNodesToDraw.clear();
	    		for(CNetworkNode currentNode : nodesToDraw){
	    			nextNodesToDraw.addAll(currentNode.DrawEdges(noOfTaxa));
	    		}
	    		//clear the array ready for the next lot
	    		nodesToDraw.clear();
	    	}
	    	while(nextNodesToDraw.size()>0);
    	}
    }
}