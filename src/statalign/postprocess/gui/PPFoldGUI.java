package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;

import javax.swing.JPanel;
import javax.swing.JScrollPane;

public class PPFoldGUI extends JPanel {

	/**
	 * This class implements the graphical interface for showing the base-pairing matrix of the current 
	 * consensus structure generated by PPFold.
	 * 
	 * @author Preeti Arunapuram
	 */
	private static final long serialVersionUID = 1L;

	public String title;
	
	private int dimension;
	private int x;
	private int y;
	private Graphics2D gr2;
	private JScrollPane scroll;
	
	public static final int OFFSET = 5;
	private int count = 0;
	
	static final int COLUMN_WIDTH = 9;
	static final int FONT_HEIGHT = 15;
	static final int OFFSET_X = 10;
	static final int OFFSET_Y = 10;
	static final int TITLE_Y = 20;
	
	
	private ColorGradient cg = new ColorGradient(Color.WHITE, Color.BLUE);
	
	private float[][] probs;
	
	
	public PPFoldGUI(String title, JScrollPane scroll) {
		this.title = title;
		this.scroll = scroll;
		dimension = 0;

		x = 10*OFFSET_X;
		y = 5*OFFSET_Y;
		
		this.setBackground(Color.BLACK);
	}
	
	protected void paintComponent(Graphics gr) {
		super.paintComponent(gr);
		gr2 = (Graphics2D)gr;
		gr2.setColor(Color.WHITE);
		if(probs == null) {
			gr2.drawString("Waiting for data..", 30, 30);
            return;
		}
		
		gr2.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		gr2.drawString(title, OFFSET_X, TITLE_Y);
		gr2.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		gr2.clearRect(x, y, dimension, dimension);
		
		
		for(int i = x; i < x+dimension; i = i+OFFSET) {
			for(int j = y; j < y+dimension; j = j+OFFSET) {
				float probability = probs[(i-x)/OFFSET][(j-y)/OFFSET];
				probability = 1/(float)Math.abs(Math.log(probability));
				gr2.setPaint(cg.getColor(probability));
				gr2.fillRect(i, j, OFFSET, OFFSET);
			}
		}
		
		count++;
		
	}
	
	public void setMatrix(float[][] probMatrix) {
		probs = probMatrix;
	}
	
	public int getDimension() {
		return dimension;
	}
	
	public void changeDimension(int x) {
		dimension = x;
	}
	
	/**
	 * This function tells the minimum size of the panel
	 */
	public Dimension getMinimumSize(){
		if(this.getPreferredSize().getHeight() > this.getHeight()) {
			scroll.createVerticalScrollBar();
		}
		
		if(this.getPreferredSize().getWidth() > this.getWidth()) {
			scroll.createHorizontalScrollBar();
		}
		
		return getPreferredSize();
	}
	
	public Dimension getPreferredSize() {
		int xPreferredSize = 20*OFFSET_X + dimension + 30;
		int yPreferredSize = 2*TITLE_Y + dimension + 30;
		
		if(yPreferredSize > this.getHeight()) {
			System.out.println("Vertical bar should be working!");
			scroll.createVerticalScrollBar();
		}
		
		if(xPreferredSize > this.getWidth()) {
			scroll.createHorizontalScrollBar();
			System.out.println("Horizontal bar should be working");
		}
		
		return new Dimension(20*OFFSET_X + dimension + 30, 2*TITLE_Y + dimension + 30);
//		if (alignment != null && alignment[0] != null) {
//			return new Dimension((alignment[0].length() + 3) * COLUMN_WIDTH + 6 * OFFSET_X, 100);
//		} else {
//			return new Dimension(0,100);
//		}
	}
	
}
	