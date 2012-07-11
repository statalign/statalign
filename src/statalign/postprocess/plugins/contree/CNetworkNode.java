package statalign.postprocess.plugins.contree;

import java.util.ArrayList;
/**
 * Nodes in the consensus network
 * Also includes a substantial amount of code that draws the edges from the node.
 * 
 * @author wood
 *
 */
public class CNetworkNode {
	// List of edges that join the node to other nodes
	public ArrayList<CNetworkEdge> joins;
	// if it is a taxon... need a name.
	public String Taxaname = null;
	// if printing the output, stores a number for the nexus file format
	public int outputNumber = 0;
	// node that this node was copied to when adding a certain split
	public CNetworkNode CopiedToNode;
	// paths this node is involved in during adding a split.  Used to find other paths to link up to.
	public ArrayList<CNetworkPath> paths;
	// flag for if considered in a function to find the place to add a compatible split
	public boolean consideration;
	// flag to show if drawn or not (i.e. assigned a position)
	public boolean drawn;
	// direction of centre of angle available to add edges in
	public double angleDirection;
	// width of angle available to add edges in
	public double angleWidth;
	// position of node in graph
	public double xPos,yPos;
    /**
     * Initialise the node
     */
	public CNetworkNode(){
		joins = new ArrayList<CNetworkEdge>();
		paths = new ArrayList<CNetworkPath>();
		consideration = false;
		CopiedToNode = null;
		outputNumber = 0;
		drawn = false;
		angleWidth = 0.0;
		// for testing...
		xPos = -1.0;
		yPos = -1.0;
	}
    /**
     * Add a path to the node path list..
     * 
     * @param pathToAdd Path to add to this node's list
     * @param nodesToClearUp List of nodes to remove paths from at end of adding split
     */
	public void AddToPaths(CNetworkPath pathToAdd,ArrayList<CNetworkNode> nodesToClearUp){
		if(paths.isEmpty()==true){
			nodesToClearUp.add(this);
		}
		paths.add(pathToAdd);
	}
    /**
     * Recursive function used to find node at which we can add a split compatible with the network
     * Goes up from first node and finds where edges joining are no longer subset of specified side of split we want to add
     * 
     * @param addSplit Split that we are adding
     * @param consideredNodes Nodes that we have so far considered
     * @param zeroSide True if looking on zero side of split, otherwise false.  Try to choose side with fewer taxa to improve efficiency. 
     * @param noOfTaxa Number of taxa
     */
	public CNetworkNode findFirstNonSubset(CNetworkSplit addSplit, ArrayList<CNetworkNode> consideredNodes,boolean zeroSide, int noOfTaxa){
	// new method to check if is subset or not...
		// note down that we are now using this node...
		consideredNodes.add(this);
		consideration = true;
		// Look at all joining edges
		for(CNetworkEdge currentEdge :joins){
			// keep looking if still a subset...
			if(addSplit.isSubset(currentEdge.split,zeroSide,noOfTaxa)==true){
				//split is not final one, then continue and find the next targets....
				// Temp array to store the joining nodes
				CNetworkNode[] joiningEdges = new CNetworkNode[]{currentEdge.networkNodeA,currentEdge.networkNodeB};
				for(CNetworkNode currentNode: joiningEdges){
					// Only go to a joining node if we haven't been there before
					if(currentNode.consideration==false){
						//Keep going...
						CNetworkNode testNode = currentNode.findFirstNonSubset(addSplit,consideredNodes,zeroSide, noOfTaxa);
						if(testNode!=null){
							return testNode;
						}
					}
				}
			}
			else{
				for(CNetworkNode consideredNode : consideredNodes){
					consideredNode.consideration = false;
				}
				consideredNodes.clear();
				return this;
			}		
		}
		//return a null if we have included all the nodes joining this (i.e. it is a taxon and maybe other cases...)
		return null;
	}
    /**
     * Function to draw the network i.e. find x and y co-ordinates for nodes.
     *  
     * @param noOfTaxa Number of taxa
     */
	public ArrayList<CNetworkNode> DrawEdges(int noOfTaxa){
		// Array list to store the nodes that we return...
		ArrayList<CNetworkNode> nodesToReturn = new ArrayList<CNetworkNode>();
		// create a 2D arraylist of the joins...
		ArrayList<ArrayList<CNetworkEdge>> orderedJoins = new ArrayList<ArrayList<CNetworkEdge>>();
		// first get the weightings and the order...
		double totalWeighting = 0.0;
		// create a list of the drawn edges in the order they occur...
		// THIS ORDER WILL BE REVERSED TO THEIR DIRECTIONS GIVEN!!
		orderedJoins.add(new ArrayList<CNetworkEdge>());
		orderedJoins.add(new ArrayList<CNetworkEdge>());
		for(CNetworkEdge currentEdge : joins){
    		// if already drawn we want to add it to a list that we can add the others to...
    		if(currentEdge.drawn == true || currentEdge.split.direction >=0.0){
    			// form first entry in box list based on the directions that have been stored in the splits
    			// may not contain anything...and that should be fine!
    			// Could  be optimised as a do while...?
    			// list 0 is for those less than this node's angle direction, list 1 is for those above it.
    			double currentDirection = (currentEdge.drawn==true)?returnAngle(currentEdge.split.direction+Math.PI):currentEdge.split.direction;
    			int j = (currentDirection<angleDirection)?0:1;
    			int storedPos = 0;
    			for(int i=0;i<orderedJoins.get(j).size();i++){
    				if(currentDirection>returnAngle(orderedJoins.get(j).get(i).split.direction+((orderedJoins.get(j).get(i).drawn==true)?Math.PI:0.0))){
    					storedPos++;
    				}
    			}
    			orderedJoins.get(j).add(storedPos,currentEdge);
    			currentEdge.considered = true;
    		}
    	}
    	orderedJoins.get(1).addAll(orderedJoins.get(0));
    	orderedJoins.get(0).clear();
   		// ONLY draw an edge if it hasn't already been drawn...
    	// populate the neighbours list...
    	for(CNetworkEdge currentEdge : joins){
    		if(currentEdge.drawn == false){
    			// renew list:
				currentEdge.neighbours = new ArrayList<CNetworkEdge>();
				// The zero element is the current edge
	    		// Get the other node from the edge
				CNetworkNode linkedNode = (currentEdge.networkNodeA == this)?currentEdge.networkNodeB:currentEdge.networkNodeA;
	    		// work out it's weighting by looking at splits
	    		currentEdge.weighting = linkedNode.FindWeighting(currentEdge, noOfTaxa);
	    		totalWeighting += currentEdge.weighting;
	    		// see if it is a box or not... i.e. do we need to add it to a list that already exists...?
	    		// go through all joined nodes if might need to state is a box...
	    		if(currentEdge.considered == false){
		    		for(CNetworkEdge secondOrderLink : linkedNode.joins){
		    			// if the corresponding split occurs 
		    			for(CNetworkEdge currentCompEdge : joins){
		    				if(secondOrderLink.split == currentCompEdge.split && secondOrderLink !=currentCompEdge){
		    					// If a second order edge matches a first order edge then we will soon have a box
		    					// Therefore we need to ensure that these two edges occur next to each other in the drawing of nodes
		    					currentEdge.neighbours.add(currentCompEdge);
		    		    		currentEdge.weighting = linkedNode.FindWeighting(currentEdge, noOfTaxa)/2.0;
		    		    		totalWeighting -= (currentEdge.weighting);
		    		    	}
		    			// end 2nd loop through current joins
		    			}
		    		// end loop of 2nd degree neighbours
		    		}
		    	// end if not considered 2nd degree neighbours
	    		}
	    	// end if not drawn
    		}
    	// end going through joins
    	}
    	// Let's now try and add to any list that exists already...
    	// add to the below one by looking for ends
    	// false if nothing found in this run...
    	boolean found = true;
    	// A ridiculous amount of nesting....
    	if(orderedJoins.get(1).size()>0){
	    	do{
	    		found = false;
	    		// probably speeded up a little bit with Do loop...?!
		    	for(CNetworkEdge currentEdge : joins){
		    		if(currentEdge.considered == false && currentEdge.drawn == false){
		    			for(CNetworkEdge neighbourLink : currentEdge.neighbours){
		    				// case is higher end, then add it to this side
		    				if(neighbourLink == orderedJoins.get(1).get(0)){
		    					orderedJoins.get(1).add(0,currentEdge);
		    					currentEdge.considered = true;
		    					found = true;
		    				}
		    				// case is lower end, then add it to position "zero"
		    				else if(neighbourLink == orderedJoins.get(1).get(orderedJoins.get(1).size()-1)){
		    					orderedJoins.get(1).add(currentEdge);
		    					currentEdge.considered = true;
		    					found = true;
		    				}
		    			// end for loop through neighbours to the join we are considering (ARGH!)
		    			}
		    		// end checking is not drawn or added as a neighbour yet...
		    		}
		    	// end loop through joins
		    	}
		    // end DO WHILE...
	    	}
	    	while(found == true);
	    // end condition of only bothering to add to the drawn nodes if nodes are already drawn...
    	}
    	// Now let's go through the rest and link them up...
    	for(CNetworkEdge neighbourLink:joins){
    		if(neighbourLink.considered == false && neighbourLink.drawn == false){
    			ArrayList<CNetworkEdge> buildList = new ArrayList<CNetworkEdge>();
    			neighbourLink.considered = true;
    			// add on start side...
    			ListToAdd(neighbourLink,true,buildList);
    			// add this node...
    			buildList.add(neighbourLink);
    			// add on other end...
    			ListToAdd(neighbourLink,false,buildList);
    			orderedJoins.add(buildList);
    		}
    	}
    	// now let's create ordered lists of the remaining ones from the neighbour array...
    	
    	// remove all the considered flags..
    	for(CNetworkEdge currentEdge : joins){
    		currentEdge.considered = false;
    	}
    	
    	
    	
    	
    	// Now we have all the joins then let's divi up the space we have and assign directions...
    	double weightMultiplier = angleWidth/totalWeighting;
    	// Now we rearrange the lists so that the final one is for all above the direction and "1" is for all below...
    	// Temporary arrays to store above and below directed assigned area.
    	ArrayList<CNetworkEdge> lessThanEdges = new ArrayList<CNetworkEdge>();
    	ArrayList<CNetworkEdge> moreThanEdges = new ArrayList<CNetworkEdge>();
    	// REMINDER: list 0 is for those less than direction, list 1 for those above it.
        boolean found2 = false;
    	for(int k=0;k<orderedJoins.get(1).size();k++){
    		if(orderedJoins.get(1).get(k).drawn == true){
    			found2 = true;
    		}
    		else{
    			if(found2 == true){
    				lessThanEdges.add(orderedJoins.get(1).get(k));
        		}
    			else{
    				moreThanEdges.add(orderedJoins.get(1).get(k));
    			}
    		}	
    	}
        // now position these correctly in the orderedlist...
        orderedJoins.add(moreThanEdges);
        orderedJoins.get(1).clear();
        orderedJoins.set(1,lessThanEdges);
    	// Need to consider case there is nothing in the first list (first node)
    	// once again go through the edges but this time in order of the lists that we created...
        // assign the beginning angle:
        double currentMinAngle = returnAngle(angleDirection - angleWidth/2);
        for(int i=0;i<orderedJoins.size();i++){
        	for(int j=0;j<orderedJoins.get(i).size();j++){
        		// no we add in order
        		CNetworkEdge currentJoin = orderedJoins.get(i).get(j);
        		// only draw those that need to be drawn
        		if(currentJoin.drawn == false){
        			// case does not have direction means we need to calculate direction...
    				double widthToAdd = weightMultiplier*currentJoin.weighting;
        			double directionToUse = returnAngle(currentMinAngle+widthToAdd/2.0);
        			if(currentJoin.split.direction<0){
        				currentJoin.split.direction = directionToUse;
        				// the following hackery is used to make sure it doesn't draw boxes that are too big...
        				double maxAngle = 2.5;
        				// if it has a neighbour, make sure it is within the limits
        				if(j>0){
        					// second angle is ALWAYS further round in +ve direction than latter....
        					// case new angle is bigger than last.... simply subtract to find difference
        						double boxAngleWidth = currentJoin.split.direction - orderedJoins.get(i).get(j-1).split.direction;
        					// otherwise difference is 2PI subtract this...
        					if(currentJoin.split.direction < orderedJoins.get(i).get(j-1).split.direction){
        						boxAngleWidth = 2*Math.PI-boxAngleWidth;
        					}
        					// now if this is > 3 then limit it...
        					if(boxAngleWidth > maxAngle){
        						currentJoin.split.direction = returnAngle(orderedJoins.get(i).get(j-1).split.direction+maxAngle);
        					}
        				}
        				// maybe there is a neighbour above...
        				else if(j==0 && orderedJoins.get(i).size() > 1){
        					if(orderedJoins.get(i).get(1).split.direction>=0.0){
        						// next angle is ALWAYS further round in +ve direction than latter....
            					// case new angle is bigger than last.... simply subtract to find difference
            						double boxAngleWidth = orderedJoins.get(i).get(1).split.direction - currentJoin.split.direction;
            					// otherwise difference is 2PI subtract this...
            					if(currentJoin.split.direction > orderedJoins.get(i).get(1).split.direction){
            						boxAngleWidth = 2*Math.PI-boxAngleWidth;
            					}
            					// now if this is > 3 then limit it...
            					if(boxAngleWidth > maxAngle){
            						currentJoin.split.direction = returnAngle(orderedJoins.get(i).get(1).split.direction-maxAngle);
            					}
        					}
        					
        				}

        			}
        			// Now draw in the edge & node given the direction and assign the correct amount of angle to the node.
        			//Draw the edge...
        			currentJoin.xPosA = xPos;
        			currentJoin.yPosA = yPos;
        			currentJoin.xPosB = xPos+currentJoin.split.edgelength*Math.sin(currentJoin.split.direction);
        			currentJoin.yPosB = yPos+currentJoin.split.edgelength*Math.cos(currentJoin.split.direction);
        			currentJoin.drawn = true;
        			// Draw the node...
        			CNetworkNode linkedNode = (currentJoin.networkNodeA == this)?currentJoin.networkNodeB:currentJoin.networkNodeA;
        			linkedNode.xPos = currentJoin.xPosB;
        			linkedNode.yPos = currentJoin.yPosB;

        			if (linkedNode.drawn == true){
        				//find smallest angle beginning, should be mutually exclusive...:
        				double drawnMax = returnAngle(linkedNode.angleDirection+(linkedNode.angleWidth/2.0));
        				double drawnMin = returnAngle(linkedNode.angleDirection-(linkedNode.angleWidth/2.0));
        				double toAddMax = returnAngle(directionToUse+(widthToAdd/2.0));
        				double toAddMin = returnAngle(directionToUse-(widthToAdd/2.0));
        				if(Math.pow((drawnMax-toAddMin),2)<0.0001){
        					linkedNode.angleDirection = returnAngle(linkedNode.angleDirection+(widthToAdd/2.0));
        					//add half the mean angle width to the smaller one...
        				}
        				else{
        					linkedNode.angleDirection = returnAngle(directionToUse+(linkedNode.angleWidth/2.0));
        				}
        			}
        			else{
        				linkedNode.angleDirection = directionToUse;
        			}
    				linkedNode.angleWidth += widthToAdd;
        			/*System.out.println("*** This Node:"+linkedNode.Taxaname);
        			System.out.println("*** Direction:"+linkedNode.angleDirection);
        			System.out.println("*** Width:    "+linkedNode.angleWidth);
        			System.out.println("*** Weight:    "+currentJoin.weighting);
        			System.out.println("*** Total Weight:    "+totalWeighting);
        			*/
        			linkedNode.drawn = true;
        			// make sure that we run the algorithm on this too!
        			nodesToReturn.add(linkedNode);
        			// Add to the current min angle...
        			currentMinAngle += widthToAdd;
        		}
        	}
        	
        }
        return nodesToReturn;
	}
	
    /**
     * Finds the weighting of this this node using the joins that join to it.  Based on:
     * a) returning "1" if it is a taxon.
     * b) return the number of linked taxon in this direction by finding which node is a subset of this split
	 * 
     * @param compEdge Edge that we want to find the weighting for
     * @param noOfTaxa Number of taxa
     */
	// Function to work out the weighting of this node by:
	public double FindWeighting(CNetworkEdge compEdge,int noOfTaxa){
		// case is a taxon
		if (joins.size() ==1) return 1;
		else{
			// go through all the edges and look for one that hasn't been drawn...
			for(CNetworkEdge currentEdge : joins){
				// Don't consider edges that have already been drawn...
				if(currentEdge.drawn == false && currentEdge != compEdge){
					// So for this edge that hasn't been drawn yet, see if it is compatible or not...
					if(currentEdge.split.isCompatible(compEdge.split,noOfTaxa)== true){
						//Find the side that one is a subset of and return the number with this side as the weighting.
						return compEdge.split.getSubsetJoins(currentEdge,noOfTaxa);
					}
				}
			}
			// if still haven't returned anything, then everything has been drawn and so return something very small (otherwise later division by zero)!...
			return 0.05;
		}
	}
	
    /**
     * Function returns joined up ArrayList if must be joined in a certain order...
     * Used for stitching together neighbours
     * Runs recursively!
     * 
     * @param neighbourLink 
     * @param start begin from start if true, begin from end if false
     * @param ListToAddTo
     */
	// Join across from start or end...
	public void ListToAdd(CNetworkEdge neighbourLink, boolean start,ArrayList<CNetworkEdge> ListToAddTo){
		// create a new list to add to and return...
		// only work on the nodes we can...
		for(CNetworkEdge neighbour : neighbourLink.neighbours){
			// if not drawn yet (shouldn't be anyway...) and not considered (not previous in chain)... then add
				if(neighbour.considered == false && neighbour.drawn == false){
					// Add to the list in the correct place
					ListToAddTo.add((start==true)?0:ListToAddTo.size()-1,neighbour);
					// mark as considered
					neighbour.considered = true;
					// now look next to it...
					ListToAdd(neighbour, start, ListToAddTo);
					// and head straight back up ignoring any subsequent neighbours as they shouldn't exist (hopefully)...
					return;
			}
		}
	// No neighbours.......... 
		return;
	}
    /**
     * Return an angle that is in the range 0 to 2*pi....
     *  
     * @param inputAngle Any input angle
     */
	public double returnAngle(double inputAngle){
		inputAngle %= 2*Math.PI;
		if(inputAngle < 0)
			inputAngle += 2*Math.PI;
		return inputAngle;
	}
	

}