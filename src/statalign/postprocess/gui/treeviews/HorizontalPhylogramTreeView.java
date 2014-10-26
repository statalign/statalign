package statalign.postprocess.gui.treeviews;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.ImageIcon;
import javax.swing.JToggleButton;

import statalign.postprocess.plugins.TreeNode;

/**
 * The graphical interface for showing the current tree.
 * @author Jonsson based on previous work by miklos, novak
 */
public class HorizontalPhylogramTreeView extends TreeView {

    // Variables

    private int noOfSamples;

    // Functions

    @Override
    public JToggleButton getToolBarButton() {
        JToggleButton button = new JToggleButton(new ImageIcon(
        		ClassLoader.getSystemResource("icons/phylogram.png")));
        button.setToolTipText("Tree visualized with branch lengths");
        return button;
    }

    @Override
    public void newSample(TreeNode root) {
        super.newSample(root);
        noOfSamples++;
    }

    @Override
    public void beforeFirstSample(int noOfTaxa) {
        super.beforeFirstSample(noOfTaxa);
        root = null;
        this.repaint();
        noOfSamples = 1;
    }

    /** It repaints the graphics */
    @Override
    public void paintComponent(Graphics gr) {
        if (root == null) {
        	return;
        }
    	
        //System.out.println("Updating tree");
        //String newick = owner.mcmc.tree.printedTree();
        Graphics2D g = (Graphics2D) gr.create();

        g.setBackground(Color.WHITE);

        g.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
        //g.setClip(0,0, panel.getWidth(), panel.getHeight());

        g.clearRect(0, 0, getWidth(), getHeight());
        g.setColor(Color.BLACK);
        //g.drawString(newick,100,100);
        //CTreeNode root = owner.output.getRoot();

        if (root == null) {
            g.drawString("Waiting for 2 sample trees...", 30, 30);
            return;
        }

        root.countLeaves();
        double maxDepth = root.maxDepth() * 1.15;

        drawHorzTree(g,
                root,
                maxDepth,
                (int) (0 + getHeight() * 0.05), (int) (getHeight() * 0.95),
                (int) (getWidth() * 0.05));
    }

    private void drawHorzTree(Graphics g, TreeNode v, double maxDepth, int y1, int y2, int x) {
        if (v.isLeaf()) {
            String s = ((v.name.trim()).length() > 20 ? v.name.substring(0, 15) + "..." : v.name);
            g.drawString(s, x + 5, (y1 + y2) / 2);
            g.drawLine(x - (int) (v.edgeLength * getWidth() * 0.9 / maxDepth), (y1 + y2) / 2, x, (y1 + y2) / 2);
        } else {
            int nrChildren = v.children.size();

            //divide the available vertical space between each child based on leafs to account for
            int[] ySplits = new int[nrChildren];
            ySplits[nrChildren - 1] = y2; //use y2 as the bottom split

            for (int i = 0; i < nrChildren - 1; i++) {
                if (i == 0) {
                    ySplits[i] = y1 + (y2 - y1) * v.children.get(i).leafCount / (v.leafCount);
                } else {
                    ySplits[i] = ySplits[i - 1] + (y2 - y1) * v.children.get(i).leafCount / (v.leafCount);
                }
            }

            if (v.edgeLength > 0.0d) {
                g.drawLine(x, ySplits[0], x - (int) (v.edgeLength * getWidth() * 0.9 / maxDepth), ySplits[0]);
                if (v.hasProperty("noOfOccurrences")) {
                    Integer noOfOccurrences = (Integer) v.getProperty("noOfOccurrences");
                    String s = String.format("%.0f", noOfOccurrences * 100 / (double) noOfSamples);
                    g.drawString(s, x + 5, ySplits[0] + 4);
                }
            }

            //place branch positions for the edge leading to each child
            int[] yBranches = new int[nrChildren];
            for (int i = 0; i < nrChildren; i++) {
                if (i == 0) {
                    yBranches[i] = (!v.children.get(i).isLeaf()) ?
                            y1 + (ySplits[i] - y1) * v.children.get(i).children.get(0).leafCount / (v.children.get(i).leafCount) :
                            y1 + (ySplits[i] - y1) / 2;

                } else {
                    yBranches[i] = (!v.children.get(i).isLeaf()) ?
                            ySplits[i - 1] + (ySplits[i] - ySplits[i - 1]) * v.children.get(i).children.get(0).leafCount / (v.children.get(i).leafCount) :
                            ySplits[i - 1] + (ySplits[i] - ySplits[i - 1]) / 2;
                }
            }
            g.drawLine(x, yBranches[0], x, yBranches[nrChildren - 1]);

            for (int i = 0; i < nrChildren; i++) {
                if (i == 0) {
                    drawHorzTree(g, v.children.get(i), maxDepth, y1, ySplits[i], x + (int) (v.children.get(i).edgeLength * getWidth() * 0.9 / maxDepth));

                } else {
                    drawHorzTree(g, v.children.get(i), maxDepth, ySplits[i - 1], ySplits[i], x + (int) (v.children.get(i).edgeLength * getWidth() * 0.9 / maxDepth));
                }
            }
        }
    }

    /** Gives the minimum size of the component */
    @Override
    public Dimension getMinimumSize() {
        return getPreferredSize();
    }

    /** It gives the preferred size of the component */
    @Override
    public Dimension getPreferredSize() {
        return new Dimension(0, noOfTaxa * 30);
    }

}
