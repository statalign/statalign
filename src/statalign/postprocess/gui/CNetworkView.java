package statalign.postprocess.gui;

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

import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.postprocess.plugins.contree.CNetwork;
import statalign.postprocess.plugins.contree.CNetworkEdge;
import statalign.postprocess.plugins.contree.CNetworkNode;
import statalign.postprocess.plugins.contree.CNetworkSplit;
import statalign.postprocess.utils.DisplayString;

/**
 * The graphical interface for showing the current network.
 *
 * @author wood, novak
 *
 */
public class CNetworkView extends JPanel implements MouseWheelListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	// variables for using in storing network data to enable resizing and zooming
	public ArrayList<DisplayString> leaflabel;
	public ArrayList<Integer[]> lines;
	// network extent in construction space
	public int minX,minY,maxX,maxY;
	// position of visible screen on intial zoom
	public int visX,visY;
	// zoom  factor
	public int zoomFactor=20;
	// Graphic
	public Graphics2D g;
	// Boolean if stored directions in array or not (saves recalculating twice)
	public boolean putInArray = false;
	// Graphic window variables
	public int height;
	public int width;
	public JScrollPane scroll;
	public final int zoomDivisor = 20;
	// Input network!
	public CNetwork network;

    // Added to comply to the TreeView interface
    /** 
     * Initialise network viewer
     * 
     * @param scroll Pane into which to place network graphic
     * 
     * */
    public CNetworkView(JScrollPane scroll) {
        this.addMouseWheelListener(this);
        this.scroll = scroll;
        lines = new ArrayList<Integer[]>();
        leaflabel = new ArrayList<DisplayString>();
    }



    // End of additions

    /**
	 * It repaints the graphics by grabbing the new data from the network.
	 * 
	 * @param gr Graphic to repaint
	 */
	@Override
	public void paintComponent(Graphics gr){
		// clear arrays
		leaflabel.clear();
		lines.clear();
		g = (Graphics2D)gr.create();
		g.setBackground(Color.WHITE);
	    g.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		g.clearRect(0, 0, this.getWidth(), this.getHeight());
	    g.setColor(Color.BLACK);
	    // if no data yet...
        if (network == null) {
            g.drawString("Waiting for data..", 30, 30);
            return;
        }
        if (putInArray == false){
			// 2D array to store the lines we want to draw
			// 2nd dimension: 0=x1,1=y1,2=x2,3=y2
			//TODO: is it more efficient to do this as 4 separate lists?
        	// Add in the edges
			for(CNetworkSplit currentSplit : network.splits){
				for(CNetworkEdge currentEdge : currentSplit.edges){
					lines.add(new Integer[]{(int)(1000*currentEdge.xPosA),(int)(1000*currentEdge.yPosA),(int)(1000*currentEdge.xPosB),(int)(1000*currentEdge.yPosB)});
				}
			}
			// Make sure we scale it to fit screen by setting initial max & min values
			if(network.nodes.size()>0){
				minX=(int)(1000*network.nodes.get(0).xPos);
				minY=(int)(1000*network.nodes.get(0).yPos);
				maxX=(int)(1000*network.nodes.get(0).xPos);
				maxY=(int)(1000*network.nodes.get(0).yPos);
			}
			// Add in the nodes
			for(CNetworkNode currentNode : network.nodes){
				// update max & min values
				minX = ((int)(1000*currentNode.xPos)<minX)?(int)(1000*currentNode.xPos):minX;
				maxX = ((int)(1000*currentNode.xPos)>maxX)?(int)(1000*currentNode.xPos):maxX;
				minY = ((int)(1000*currentNode.yPos)<minY)?(int)(1000*currentNode.yPos):minY;
				maxY = ((int)(1000*currentNode.yPos)>maxY)?(int)(1000*currentNode.yPos):maxY;
				// Remove if statement below to leave in internal node names...
				if(currentNode.joins.size()==1){
					String printName = (currentNode.Taxaname.length()> 10)?currentNode.Taxaname.substring(0, 10)+"...":currentNode.Taxaname;
					DisplayString addNow = new DisplayString(printName,(int)(1000*currentNode.xPos),(int)(1000*currentNode.yPos));
					leaflabel.add(addNow);
				}
			}
        }
        redraw();
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
		height = (this.getHeight()-30);
		width = (this.getWidth()-100);
		double xMultiple =  (double)Math.abs((double)width/(double)(maxX-minX));
		double yMultiple =  (double)Math.abs((double)height/(double)(maxY-minY));
		int xBorder = 10;
		int yBorder = 20;
		for(int i=0;i<lines.size();i++){
			Integer[] curLine = lines.get(i);
			g.drawLine((int)((curLine[0]-minX)*xMultiple)+xBorder,(int)((curLine[1]-minY)*yMultiple)+yBorder,(int)((curLine[2]-minX)*xMultiple)+xBorder,(int)((curLine[3]-minY)*yMultiple)+yBorder);
		}
		for(DisplayString curString:leaflabel){
			g.drawString(curString.label,(int)((curString.x-minX)*xMultiple)+xBorder,(int)((curString.y-minY)*yMultiple)+yBorder);
		}
	}

	/**
	 * It listens for mousewheel movement and sets the drawing appropriately.
	 */
    public void mouseWheelMoved(MouseWheelEvent e) {
        int oldZoomFactor = zoomFactor;
        zoomFactor -= e.getWheelRotation();
        zoomFactor = (zoomFactor > zoomDivisor * 3) ? zoomDivisor * 3 : zoomFactor;
        zoomFactor = (zoomFactor < zoomDivisor) ? zoomDivisor : zoomFactor;
        // calculate all the positions in the unscaled co-ords:
        double xMousePos = ((double) (scroll.getViewport().getViewPosition().x + MouseInfo.getPointerInfo().getLocation().x - scroll.getLocationOnScreen().x)) / ((double) oldZoomFactor / (double) zoomDivisor);
        double yMousePos = ((double) (scroll.getViewport().getViewPosition().y + MouseInfo.getPointerInfo().getLocation().y - scroll.getLocationOnScreen().y)) / ((double) oldZoomFactor / (double) zoomDivisor);
        double xNewView = (xMousePos - (double) scroll.getWidth() / (2 * (double) zoomFactor / (double) zoomDivisor));
        double yNewView = (yMousePos - (double) scroll.getHeight() / (2 * (double) zoomFactor / (double) zoomDivisor));
        xNewView = (xNewView < 0) ? 0 : xNewView;
        yNewView = (yNewView < 0) ? 0 : yNewView;
        xNewView = (xNewView > scroll.getWidth()) ? scroll.getWidth() - (scroll.getWidth() / zoomFactor * zoomDivisor) : xNewView;
        yNewView = (yNewView > scroll.getHeight()) ? scroll.getHeight() - (scroll.getHeight() / zoomFactor * zoomDivisor) : yNewView;
        scroll.getViewport().setViewPosition(new Point((int) (xNewView * zoomFactor / (double) zoomDivisor), (int) (yNewView * zoomFactor / (double) zoomDivisor)));
        scroll.repaint();
        scroll.setViewportView(this);
    }

}
