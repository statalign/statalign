package statalign.postprocess.plugins;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import statalign.base.Vertex;

public class TreeNode {

    // Variables

    /** The name/identifier of this node. */
    public String name;

    /** The branch length of this node. */
    public double edgeLength;

    /** The children of this node. */
    public List<TreeNode> children;

    /** The parent of this node. */
    public TreeNode parent;

    /** Properties associated with this node. Lazy initialized: Not instantiated unless needed. */
    public HashMap<String, Object> properties;
    
    // TODO: switch out for a property.
    public boolean ignoreNode = false;

    // Functions

    /**
     * A constructor which does not specify a name nor a edge length.
     */
    public TreeNode() {
        this(null, 0.0d);
    }
    
//    /**
//     * Copy constructor
//     */
//    public TreeNode(TreeNode other) {
//        this(other.name, other.edgeLength);
//        System.out.print(this+"="+name);
//        if (other.isLeaf()) {
//        	children.clear();
//        }
//        else {          	
//	        for (int i=0; i<children.size(); i++) {
//	        	System.out.print("(");
//	        	children.set(i,new TreeNode(other.children.get(i)));
//	        	System.out.print("):"+edgeLength);
//	        }
//        }
//    }     
    
    public TreeNode(Vertex v) {
        this(v.name, v.edgeLength);
//        System.out.print(System.identityHashCode(this)+"="+name+"["+v.leafCount+"]");           	
//        System.out.print("(");
        if (v.left != null) children.add(new TreeNode(v.left));        
        if (v.right != null) children.add(new TreeNode(v.right));
        for (int i=0; i<children.size(); i++) {
        	children.get(i).parent = this;
        }        
//        System.out.print("):"+edgeLength);	                 
    }

    /**
     * A constructor specifying the name of the node.
     * @param name the name of the node.
     */
    public TreeNode(String name) {
        this(name, 0.0d);
    }

    /**
     * A constructor specifying the name and the edge length of the node.
     * @param name the name of the node.
     * @param edgeLength the edge length of the node.
     */
    public TreeNode(String name, double edgeLength) {
        this.name = name;
        this.edgeLength = edgeLength;
        this.children = new ArrayList<TreeNode>(2);
    }

    /**
     * Determines whether this node is a leaf or not.
     * @return whether this node is a leaf or not.
     */
    public boolean isLeaf() {
        return children.size() == 0;
    }

    /**
     * Adds a child to this tree's list of children.
     * @param node the child to be added.
     */
    public void addChild(TreeNode node) {
        children.add(node);
    }

    /**
     * Finds all of the leaves of the subtree of this node.
     * @return a list of the leaves of this node's subtree.
     */
    public List<TreeNode> getLeaves() {
        return getLeaves(this);
    }

    /**
     * Same as {@link #getLeaves()}.
     * @param node the node that we want to determine the leaves of.
     * @return a list of the leaves of this nodes subtree.
     */
    public List<TreeNode> getLeaves(TreeNode node) {
        List<TreeNode> leaves = new ArrayList<TreeNode>();

        if (node.isLeaf()) {
            leaves.add(node);
        }
        for (TreeNode n : node.children) {
            leaves.addAll(getLeaves(n));
        }

        return leaves;
    }

    /**
     * Constructs a Newick string from this subtree.
     * @return a string representation of this subtree in a Newick format.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        if (!children.isEmpty()) {
            sb.append('(').append(children.get(0));
            for (int i = 1; i < children.size(); i++)
                sb.append(',').append(children.get(i));
            sb.append(')');
        } else {
        	sb.append(name);
        }
        if (parent != null) {
            sb.append(":").append(String.format(Locale.US, "%.5f", edgeLength));
        } else {
        	sb.append(";");
        }

        return sb.toString();
    }
    
    /**
     * Constructs a Newick string from this subtree.
     * @return a string representation of this subtree in a Newick format.
     */
    public String toStringWithProbs(int noOfSamples) {
        StringBuilder sb = new StringBuilder();
        if (!children.isEmpty()) {
            sb.append('(').append(children.get(0).toStringWithProbs(noOfSamples));
            for (int i = 1; i < children.size(); i++)
                sb.append(',').append(children.get(i).toStringWithProbs(noOfSamples));
            sb.append(')');
            String s = "";
            if (hasProperty("noOfOccurrences")) {                
            	Integer noOfOccurrences = (Integer) getProperty("noOfOccurrences");            	
             	s = String.format("%.0f", noOfOccurrences * 100 / (double) noOfSamples);
            }                                   
            sb.append(s);
        } else {
        	sb.append(name);
        }
        if (parent != null) {
            sb.append(":").append(String.format(Locale.US, "%.5f", edgeLength));
        } else {
        	sb.append(";");
        }

        return sb.toString();
    }


    // Properties functions

    /**
     * Adds a key=value property to be associated with this node.
     * @param key the key portion of the property.
     * @param value the value portion of the property.
     */
    public void addProperty(String key, Object value) {
        // Lazy initialization
        if (properties == null) {
            properties = new HashMap<String, Object>();
        }
        properties.put(key, value);
    }

    /**
     * Return the value for a certain key property associated with this node.
     * @param key the key portion of the property.
     * @return the value associated with this key (or null if not existing).
     */
    public Object getProperty(String key) {
        if (properties == null) {
            return null;
        }
        return properties.get(key);
    }

    /**
     * Same as {@link #getProperty(String)} but this method returns a value casted to an {@link Integer}.
     * @param key the key portion of the property.
     * @return the value associated with this key (casted to an {@link Integer}).
     */
    public int getIntProperty(String key) {
        return (Integer) properties.get(key);
    }

    /**
     * Determines whether a certain key exists in the properties associated with this node.
     * @param key the key of the property.
     * @return a boolean whether this key actually exists or not.
     */
    public boolean hasProperty(String key) {
        return properties != null && properties.containsKey(key);
    }

    // Misc variables/functions

    /**
     * The sum of the length of the names below this vertex. Used in the tree visualisation to decide
     * how much space this subtree needs to be able to put leaves' names nicely onto the graphical
     * interface.
     */
    public int nameLengthSum;

    /** The number of leaves in this subtree. */
    public int leafCount;

    /**
     * Calculate the nodes needed to cross to get to the most distant leaf
     * Used in tree visualisation.
     * @return steps needed to reach the most distant leaf.
     */
    public int maxSteps() {
        /*if(left == null){
              return 0;
          }
          else return 1 + Math.max(left.maxSteps(),right.maxSteps());*/

        if (isLeaf()) {
            return 0;
        }
        int max = 0;
        for (TreeNode node : children) {
            max = Math.max(max, node.maxSteps());
        }
        return 1 + max;
    }

    /**
     * This function calculates the sum of the length of the strings representing the
     * names of the leaves below this subtree. Used in tree visualisation.
     * @return the sum of the lengths of the names below this vertex.
     */
    public int countNameLengthSum() {
        return nameLengthSum = isLeaf()
                ? Math.min(10, name.indexOf(' ') == -1 ? name.length() : name.indexOf(' ')) + 2
                : children.get(0).countNameLengthSum() + children.get(1).countNameLengthSum();
    }

    /**
     * Calculated the number of leads that exists below this node
     * @return the number of leaves below the node.
     */
    public int countLeaves() {
        if (this.isLeaf()) {
            return leafCount = 1;
        }

        int sum = 0;
        for (TreeNode child : children) {
            sum += child.countLeaves();
        }
        return leafCount = sum;

    }

    /**
     * Calculate the maximum depth below on this vertex.
     * Used in tree visualisation.
     * @return maximum depth below this vertex.
     */
    public double maxDepth() {
        if (this.isLeaf()) {
            return edgeLength;
        } else {
            double max = 0;
            for (int i = 0; i < children.size(); i++) {
                if (max < children.get(i).maxDepth()) {
                    max = children.get(i).maxDepth();
                }
            }
            return edgeLength + max;
        }
    }
    
    /** TODO: WHAT? */
    public TreeNode getLeft() {
    	return children.get(0);
    }
    
    public TreeNode getRight() {
    	return children.get(1);
    }

}