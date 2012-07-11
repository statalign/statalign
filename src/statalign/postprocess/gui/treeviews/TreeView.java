package statalign.postprocess.gui.treeviews;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToggleButton;

import statalign.postprocess.plugins.TreeNode;

public abstract class TreeView extends JPanel {

    // Variables

    protected TreeNode root;
    protected int noOfTaxa;
    protected String identifier;
    protected JComponent parent;

    // Abstract functions

    public abstract JToggleButton getToolBarButton();

    // Functions

    public TreeView() {
        this.identifier = this.getClass().getSimpleName();
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    public String getIdentifier() {
        return identifier;
    }

    public void beforeFirstSample() {
    }

    public void beforeFirstSample(int noOfTaxa) {
    	this.noOfTaxa = noOfTaxa;
    }

    public void newSample(TreeNode root) {
        this.root = root;
    }

    public void setParent(JScrollPane parent) {
        this.parent = parent;
    }

}