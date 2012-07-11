package statalign.postprocess.gui.treeviews;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.MouseInfo;
import java.awt.Point;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.ArrayList;

import javax.swing.ImageIcon;
import javax.swing.JScrollPane;
import javax.swing.JToggleButton;

import statalign.postprocess.plugins.TreeNode;
import statalign.postprocess.utils.DisplayString;

/**
 * The graphical interface for showing the current network.
 *
 * @author wood, novak
 *
 */
public class NetworkTreeView extends TreeView implements MouseWheelListener {
	private static final long serialVersionUID = 1L;
	
	// variables for using in storing network data to enable resizing and zooming
	public ArrayList<DisplayString> leaflabel;
	public ArrayList<Integer[]> lines;
	// network extent in construction space
	public int minX,minY,maxX,maxY;
	// position of visible screen on initial zoom
	public int visX,visY;
	// zoom  factor
	public int zoomFactor=20;
	// Graphics info.
	public Graphics2D g;
	public boolean putInArray = false;
	public int height;
	public int width;
	public JScrollPane scroll;
	public final int zoomDivisor = 20;

    // Added to comply to the TreeView interface
    public NetworkTreeView() {
        addMouseWheelListener(this);
    }

    @Override
    public JToggleButton getToolBarButton() {
        JToggleButton button = new JToggleButton(new ImageIcon(
        		ClassLoader.getSystemResource("icons/vert.png")));
        button.setToolTipText("Network view");
        return button;
    }

    @Override
    public void newSample(TreeNode root) {
        super.newSample(root);
        putInArray = false;
    }

    @Override
    public void setParent(JScrollPane parent) {
        scroll = parent;
    }

    // End of additions

    /**
	 * It repaints the graphics
	 */
	@Override
	public void paintComponent(Graphics gr){

		if (root == null) {
			return;
		}
		// Standard repainting...
			g = (Graphics2D)gr.create();
            g.setBackground(Color.WHITE); 
            g.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		    //g.setClip(0,0, panel.getWidth()*zoomFactor, panel.getHeight()*zoomFactor);
		    
			g.clearRect(0, 0, getWidth(), getHeight());
            g.setColor(Color.BLACK); 
			if(putInArray == false){
			int nameLengthCount = root.countNameLengthSum();
			root.countLeaves();
			
			// 2D array to store the lines we want to draw
			// 2nd dimension: 0=x1,1=y1,2=x2,3=y2
			//TODO: is it more efficient to do this as 4 separate lists?
			lines = new ArrayList<Integer[]>(root.leafCount);
			// 2D array to store the strings we want to join
			// 2nd dimension 0=ref to string in owner.mcmc.tree, 1=x, 2=y
			leaflabel = new ArrayList<DisplayString>(root.leafCount);
			// Initialise max & min
			minY=0;
			minX=0;
			maxY=0;
			maxX=0;
			// calculate the positions of objects...
			drawNetwork(root,0,0,0,1);
			// actually draw the network
			redraw();
			// Note we have saved the data
			putInArray = true;
		}
		else{
			redraw();
		}
	}


	/**
	 *
	 * @param v The tree node we want to calculate the position of
	 * @param x Current x position
	 * @param y Current y position
	 * @param dir The current tree direction
	 * @param prop The proportion of 360 deg. angle available
	 */
	private void drawNetwork(TreeNode v, int x, int y, double dir, double prop){
		double curDir;
		int curX,curY;
		// Test to see if this is a node or a leaf...
		if(v.isLeaf()){
			// Prepare the string for output
			String s = ((v.name.trim()).length() > 10 ? v.name.substring(0, 10)+"..." : v.name);
			leaflabel.add(new DisplayString(s,x,y));
			// Update the maximum and minimums...
			if(x<minX) minX=x;
			else if (x>maxX) maxX=x;
			if(y<minY) minY=y;
			else if (y>maxY) maxY=y;
		}
		else{
            TreeNode left = v.children.get(0);
            TreeNode right = v.children.get(1);
			// for each vertex we want to head in this direction, note that y=0 is the zero direction
			curDir = dir+Math.PI*((prop/2)*((left.leafCount/v.leafCount)-1));
			// save the vertex line to our draw line list
			// times by 500 as we are storing as integers and want something that scales well
			curX = (int)(x+Math.cos(curDir)*left.edgeLength*500);
			curY = (int)(y+Math.sin(curDir)*left.edgeLength*500);
			lines.add(new Integer[]{x, y,curX ,curY});
			//check for the boundaries now to save going through the array later
			if(curX<minX) minX=curX;
			else if (curX>maxX) maxX=curX;
			if(curY<minY) minY=curY;
			else if (curY>maxY) maxY=curY;
			drawNetwork(left,curX,curY,curDir,prop*left.leafCount/v.leafCount);
			// now do the right side of the vertex
			curDir = dir-Math.PI*((prop/2)*((right.leafCount/v.leafCount)-1));
			// save the vertex line to our draw line list
			curX = (int)(x+Math.cos(curDir)*right.edgeLength*500);
			curY = (int)(y+Math.sin(curDir)*right.edgeLength*500);
			lines.add(new Integer[]{x, y,curX ,curY});
			//check for the boundaries now to save going through the array later
			if(curX<minX) minX=curX;
			else if (curX>maxX) maxX=curX;
			if(curY<minY) minY=curY;
			else if (curY>maxY) maxY=curY;
			drawNetwork(right,curX,curY,curDir,prop*right.leafCount/v.leafCount);
		}
	}

	/**
	 * Gives the minimum size of the component
	 */

	@Override
	public Dimension getMinimumSize(){
		return getPreferredSize();
	}


	/**
	 * It gives the preferred size of the component
	 */
	@Override
	public Dimension getPreferredSize() {
		return new Dimension((scroll.getWidth()*zoomFactor/zoomDivisor)-10,(scroll.getHeight()*zoomFactor/zoomDivisor)-10);
	}
	/**
	 * It draws the network from our arrays to the correct size
	 */
	public void redraw(){
		height = (getHeight()-30);
		width = (getWidth()-100);
		double xMultiple =  Math.abs((double)width/(double)(maxX-minX));
		double yMultiple =  Math.abs((double)height/(double)(maxY-minY));
		int xBorder = 10;
		int yBorder = 20;
		for(Integer[] curLine:lines){
			g.drawLine((int)((curLine[0]-minX)*xMultiple)+xBorder,(int)((curLine[1]-minY)*yMultiple)+yBorder,(int)((curLine[2]-minX)*xMultiple)+xBorder,(int)((curLine[3]-minY)*yMultiple)+yBorder);
		}
		for(DisplayString curString:leaflabel){
			g.drawString(curString.label,(int)((curString.x-minX)*xMultiple)+xBorder,(int)((curString.y-minY)*yMultiple)+yBorder);
		}
	}


	/**
	 * Event detector for mouse wheel movement.
	 */
    public void mouseWheelMoved(MouseWheelEvent e) {
        int oldZoomFactor = zoomFactor;
        zoomFactor -= e.getWheelRotation();
        zoomFactor = (zoomFactor > zoomDivisor * 3) ? zoomDivisor * 3 : zoomFactor;
        zoomFactor = (zoomFactor < zoomDivisor) ? zoomDivisor : zoomFactor;
        //networkGUI.setBounds(0, 0, (int)((double)(networkGUI.zoomFactor/(double)networkGUI.zoomDivisor)*2*networkGUI.scroll.getWidth()), (int)((double)(networkGUI.zoomFactor/(double)networkGUI.zoomDivisor)*networkGUI.scroll.getHeight()));
        // calculate all the positions in the unscaled co-ords:
        double xMousePos = ((scroll.getViewport().getViewPosition().x + MouseInfo.getPointerInfo().getLocation().x - scroll.getLocationOnScreen().x)) / ((double) oldZoomFactor / (double) zoomDivisor);
        double yMousePos = ((scroll.getViewport().getViewPosition().y + MouseInfo.getPointerInfo().getLocation().y - scroll.getLocationOnScreen().y)) / ((double) oldZoomFactor / (double) zoomDivisor);
        double xNewView = (xMousePos - scroll.getWidth() / (2 * (double) zoomFactor / zoomDivisor));
        double yNewView = (yMousePos - scroll.getHeight() / (2 * (double) zoomFactor / zoomDivisor));
        xNewView = (xNewView < 0) ? 0 : xNewView;
        yNewView = (yNewView < 0) ? 0 : yNewView;
        xNewView = (xNewView > scroll.getWidth()) ? scroll.getWidth() - (scroll.getWidth() / zoomFactor * zoomDivisor) : xNewView;
        yNewView = (yNewView > scroll.getHeight()) ? scroll.getHeight() - (scroll.getHeight() / zoomFactor * zoomDivisor) : yNewView;
        scroll.getViewport().setViewPosition(new Point((int) (xNewView * zoomFactor / zoomDivisor), (int) (yNewView * zoomFactor / zoomDivisor)));
        scroll.repaint();
        scroll.setViewportView(this);
    }

}

