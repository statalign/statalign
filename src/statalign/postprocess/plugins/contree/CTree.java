package statalign.postprocess.plugins.contree;

import java.util.ArrayList;
import java.util.HashMap;

import statalign.postprocess.plugins.TreeNode;

public class CTree {

    // Variables

    private TreeNode root;
    public ArrayList<TreeNode> nodeList;
    public HashMap<String, Integer> parentList;
    public ArrayList<Cluster> clusters;

    // Functions

    public CTree() {
        parentList = new HashMap<String, Integer>();
        nodeList = new ArrayList<TreeNode>();

        this.root = new TreeNode("root");
        nodeList.add(this.root);
    }

    public TreeNode getRoot() {
        return root;
    }

    @Override
    public String toString() {
        return root + ";";
    }

}
