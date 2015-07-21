package statalign.postprocess.plugins.structalign;

import java.io.IOException;

import statalign.base.State;
import statalign.base.Vertex;
import statalign.model.ext.ModelExtManager;
import statalign.model.ext.ModelExtension;
import statalign.model.ext.plugins.StructAlign;
import statalign.postprocess.gui.treeviews.TreeView;
import statalign.postprocess.plugins.CurrentTree;
import statalign.postprocess.plugins.TreeNode;

public class StructTree extends CurrentTree {

	public StructAlign structAlign;
	public TreeNode node;
	
    public StructTree() {
    	screenable = true;
		outputable = true;
		postprocessable = true;
		postprocessWrite = false;
		selected = false;
		active = false;		
	}

    @Override
    public String getTabName() {
        return "Tree weighted by structural diffusivity";
    }

    @Override
    public String getTip() {
        return "Tree weighted by structural diffusivity";
    }

    @Override
    public double getTabOrder() {
        return 3.5d;
    }
    
    /**
     * Updates the tree held in this object's TreeNode.
     * This method is called from the connected 
     * StructAlign object when the tree is changed.
     */
    public void updateStructTree(Vertex root, double[] sigma2) {
    	node = new TreeNode(root);    	
//    	System.out.println();
    	updateStructTree(root,node,sigma2);
    }
    public void updateStructTree(Vertex v, TreeNode n, double[] sigma2) {
    	    	 
    	if (v.left != null) {
    		 n.children.get(0).edgeLength = sigma2[v.left.index];
    		 updateStructTree(v.left,n.children.get(0),sigma2);
    	}
    	if (v.right != null) {
    	   	n.children.get(n.children.size()-1).edgeLength = sigma2[v.right.index];
    	   	updateStructTree(v.right,n.children.get(n.children.size()-1),sigma2);
    	}
	}
    
    @Override
	public void init(ModelExtManager modelExtMan) {
    	super.init(modelExtMan);
		for(ModelExtension modExt : modelExtMan.getPluginList()) {
			if(modExt instanceof StructAlign) {
				structAlign = (StructAlign) modExt;
				structAlign.connectStructTree(this);
			}
		}
		active = structAlign.isActive() & structAlign.showStructTree;
		postprocessWrite = active;
	}
    
    @Override
    public void newPeek(State state) {
    	if(!active)
			return;
        if(show && node!=null) {        	
        	for (TreeView view : treeViews) {
        		view.newSample(node);
        	}
        	refreshGUI();
        }
    }

    @Override
    public String getFileExtension() {
        return "struct.tree";
    }

    @Override
    public void newSample(State state, int no, int total) {
    	if(!active)
			return;        
        if(show && node!=null) {
        	for (TreeView view : treeViews) {
        		view.newSample(node);
        	}

        	refreshGUI();
        }
        if (sampling) {
            try {
                if (state.isBurnin) file.write("Burnin " + no + "\tTree:\t" + node.toString() + "\n");
                else file.write("Sample " + no + "\tTree:\t" + node.toString() + "\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        if (postprocessWrite) {
            try {
                outputFile.write("[" + no + "] tree 'tree_" + no + "'= "
        				+ node.toString() + "\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

}
