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
 * @author miklos, novak
 *
 * TODO: remove?
 */
public class VerticalPhylogramTreeView extends TreeView {

	private static final long serialVersionUID = 1L;

	@Override
    public JToggleButton getToolBarButton() {
        JToggleButton button = new JToggleButton(new ImageIcon(
        		ClassLoader.getSystemResource("icons/vert.png")));
        button.setToolTipText("Vertical tree view");
        return button;
    }

    // Repaints the graphics.
    @Override
    public void paintComponent(Graphics gr) {
    	if (root == null) {
    		return;
    	}
    	
        //System.out.println("Updating tree");
        Graphics2D g = (Graphics2D) gr.create();

        g.setBackground(Color.WHITE);
        g.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
        g.setClip(0, 0, getParent().getWidth(), getParent().getHeight());
        g.clearRect(0, 0, getWidth(), getHeight());
        g.setColor(Color.BLACK);

        if (root == null) {            // TODO: move to error handler.
            g.drawString("Waiting for data..", 30, 30);
        } else {
            root.countNameLengthSum();
            double maxDepth = root.maxDepth() * 1.15;
            drawTree(g, root, maxDepth, 0, getWidth(), (int) (getHeight() * 0.05));
        }
    }

    private void drawTree(Graphics g, TreeNode node, double maxDepth, int x1, int x2, int y) {
        if (node.isLeaf()) {
            String s = ((node.name.trim()).length() > 10 ? node.name.substring(0, 10) + "..." : node.name);
            g.drawString(s, (x1 + x2) / 2 - (s.indexOf(' ') == -1 ? s.length() : s.indexOf(' ')) * 4, y + 15);
            g.drawLine((x1 + x2) / 2, y, (x1 + x2) / 2, y - (int) (node.edgeLength * getHeight() / maxDepth));
        } else {
            // HACK:
            TreeNode left = node.children.get(0);
            TreeNode right = node.children.get(1);
            int xmid = x1 + (x2 - x1) * left.nameLengthSum / (left.nameLengthSum + right.nameLengthSum);
            if (node.edgeLength > 0.0d) {
                g.drawLine(xmid, y, xmid, y - (int) (node.edgeLength * getHeight() / maxDepth));
            }
            // Best code ever
            TreeNode leftleft = left.isLeaf() ? null : left.children.get(0);
            TreeNode leftright = left.isLeaf() ? null : left.children.get(1);
            TreeNode rightleft = right.isLeaf() ? null : right.children.get(0);
            TreeNode rightright = right.isLeaf() ? null : right.children.get(1);
            int xleft = (leftleft != null) ?
                    x1 + (xmid - x1) * leftleft.nameLengthSum / (leftleft.nameLengthSum + leftright.nameLengthSum) :
                    x1 + (xmid - x1) / 2;
            int xright = (rightleft != null) ?
                    xmid + (x2 - xmid) * rightleft.nameLengthSum / (rightleft.nameLengthSum + rightright.nameLengthSum) :
                    xmid + (x2 - xmid) / 2;
            g.drawLine(xleft, y, xright, y);
            drawTree(g, left, maxDepth, x1, xmid, y + (int) (left.edgeLength * getHeight() / maxDepth));
            drawTree(g, right, maxDepth, xmid, x2, y + (int) (right.edgeLength * getHeight() / maxDepth));
        }
    }

    // Gives the minimum size of the component
    @Override
    public Dimension getMinimumSize() {
        return getPreferredSize();
    }

    // It gives the preferred size of the component
    @Override
    public Dimension getPreferredSize() {
        return new Dimension(0, 100);
    }

}