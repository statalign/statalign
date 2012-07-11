package com.ppfold.main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Stack;

import com.ppfold.algo.Node;
import com.ppfold.algo.Tree;

public class NewickReader{
	
	private static int nodeCounter = 0; //this will be used to specify node ID, 
	//which is an arbitrary number for each node mostly for debugging purposes.

	public static Tree readNewick(String fullFileName) throws Exception{
        String treestring = "";
        
        try {

			FileReader input = new FileReader(fullFileName);
			BufferedReader reader = new BufferedReader(input);
			String line;
			line = reader.readLine().trim(); //always trim leading and ending spaces off 
			if(line!=null){
				treestring = treestring.concat(line.trim()); //parse lines to one string
				line = reader.readLine();
				if(line!=null){
					treestring = treestring.concat(line.trim()); //parse lines to one string
				}
			}

			reader.close();

	        Tree tree = NewickReader.parse(treestring);
	        return tree; 
		}
		catch (IOException e){
			// If another exception is generated, print a stack trace
			throw new Exception("Error while trying to read tree file! Check that the name is OK.");
		}
		catch (Exception e){
			throw new Exception("Error while trying to parse the tree file!");
		}
	}
	
	public static Tree parse(String inputString) throws Exception {

		try{
			//remove ; from end
			if(inputString.charAt(inputString.length()-1) == ';'){
				inputString = inputString.substring(0,inputString.length()-1);
			}

			//remove any enclosing () pairs
			if(inputString.charAt(0)=='(' && (inputString.charAt(inputString.length()-1) == ')')){
				inputString = inputString.substring(1,inputString.lastIndexOf(')'));
			}

			//Create a root node
			Node root = new Node(0);
			root.setName("root");
			root.setId(0);
			root.setDistanceFromParent(0);
			nodeCounter = 1;

			//Process full string that corresponds to node=root. 
			int success = processString(inputString, root);

			if(success==0){
				//identify if there is an extra root node 
				if(root.getChildren().size()==1){
					if(root.getChildren().get(0).getDistanceFromParent() == 0d){
						root = root.getChildren().get(0); //drop own root and replace it with the one provided by user 
					}
				}
				
				Tree tree = new Tree(root); //create tree using the root node
				return tree; 
			}
			else{
				return null; 
			}
		}
		catch(Exception e){
			throw new Exception("An error occured while trying to parse the tree. Check format."); 	
		}
	
	
	}
	
	private static int processString(String inputString, Node parent) throws Exception{
		//There are two cases here. Either the inputString contains brackets (more nested nodes),
		//or it only contains comma-separated children at the last level. 
		if(inputString.indexOf('(') == -1){
			//there are no brackets at all. 
			String[] childrenString = inputString.split(","); //this will contain all the children
			
			for(String childString:childrenString){
				//for each child, extract name and distance, and add it to the parent. 
				String[] childSplitString = childString.split(":",2);
				Node thischild = new Node(Double.parseDouble(childSplitString[1])); //this sets distance from parent 
				thischild.setName(childSplitString[0]); 
				thischild.setId(nodeCounter++);
				parent.addChild(thischild);
			}
		}
		else{
			//there are brackets.
			//first identify a bracket pair, because it represents an internal node 
			//("this node" in the followings).
			//then create the node for this pair (process it recursively until needed)
			//then deal with the rest of the string. 
			int[] positionOfNodeSubstring;
			try{
				positionOfNodeSubstring = findNodeSubstring(inputString);
			}
			catch(Exception e){
				throw new Exception(e); 
			}

			
			if(positionOfNodeSubstring == null){return -1;} //error; findNodeSubstring will report what went wrong. 
			
			String nodeSubstring = inputString.substring(positionOfNodeSubstring[0]+1, positionOfNodeSubstring[1]); //the subtree coming from this node

			String leftParentString = inputString.substring(0, positionOfNodeSubstring[0]); //everything left of this node 
			String rightParentStringFull = inputString.substring(positionOfNodeSubstring[1]+1); //this contains info for this node + other nodes in parent. 
			String rightParentString = ""; //will contain everything in the parent right of this node 
			String[] splitRightParentString = rightParentStringFull.split(",",2); //isolate info about this node from other nodes in parent

			if(splitRightParentString.length==2){
				//there are extra nodes at the end
				rightParentString = splitRightParentString[1]; 
			}
			
			String nameAndNumberString = splitRightParentString[0]; //just the info about this node, always expected to be there 
			String remainderParentString = leftParentString.concat(rightParentString); //the rest of the parent
			
			//extract&set the info for this node
			String[] nameNumberSplitString = nameAndNumberString.split(":",2);			

			Node newNode = new Node();
			if(nameNumberSplitString.length==1){
				//only number is given 
				newNode.setDistanceFromParent(Double.parseDouble(nameNumberSplitString[0]));
			}
			else{
				//both name and number are given 
				newNode.setDistanceFromParent(Double.parseDouble(nameNumberSplitString[1]));
				newNode.setName(nameNumberSplitString[0]); 
			}
	
			newNode.setId(nodeCounter++);
			parent.addChild(newNode);
			
			//deal with the subtree of this node as needed
			processString(nodeSubstring, newNode);
			
			//process the rest of the parent, unless we have come to the end.
			if(remainderParentString.length()!=0){
				processString(remainderParentString, parent);
			}
			
		}
		return 0; //success

	}

	
	public static int[] findNodeSubstring(String inputString) throws Exception{
		if(inputString.indexOf('(') == -1){return null;}
        int[] pair = new int[2];
        
//		identify outermost/first pairing parantheses. 
		
		Stack<Integer> s = new Stack<Integer>();
        
        for (int i=0;i<inputString.length(); i++) {
            if (inputString.charAt(i) == '(' ) {
                s.push(i);
            }
            if (inputString.charAt(i) ==')' ) {
                if (s.isEmpty()) {
                    throw new Exception("Unmatched right paranthesis at position "+ (i+1) +" in input tree!");
                }
                int leftPos = s.pop();
                pair[0] = leftPos;
                pair[1] = i; 
            }
        }
        if (!s.isEmpty()) {
        	throw new Exception("Unmatched left paranthesis at position "+ (s.pop()+1) +" in input tree!");	
        }

		return pair; 
	}
	
	
}
