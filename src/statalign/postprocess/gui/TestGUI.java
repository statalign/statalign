package statalign.postprocess.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.MouseInfo;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.util.Random;

import javax.swing.BorderFactory;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.postprocess.plugins.PPFold;
public class TestGUI extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private JPanel panel;
	private PPFold owner;
	private JPanel plotArea;
	
	private int dimension;
	private int x;
	private int y;
	private Graphics2D gr2;
	private JScrollPane scroll;
	
	private int zoomFactor = 20;
	private int zoomDivisor = 20;
	public static final int OFFSET = 5;
	private int count = 0;
	
	
	private ColorGradient cg = new ColorGradient(Color.WHITE, Color.BLUE);
	
	private Random random = new Random();
	private float[][] probs;
	
	
	public TestGUI(PPFold owner) {
		this.owner = owner;
		dimension = 0;
		//x = (int)((0.5)*(this.getWidth()-dimension));
		x = 0;
		//y = (int)((0.5)*(this.getHeight()-dimension));
		y = 0;
		
		
		//add(scroll, BorderLayout.CENTER);
	}
	
	protected void paintComponent(Graphics gr) {
		
		
		super.paintComponent(gr);
		gr2 = (Graphics2D)gr;
		
		if(probs == null) {
			gr2.drawString("Waiting for data..", 30, 30);
            return;
		}
		
		else {
			gr2.drawString("Hello there!: " + count, 30, 30);
		}
		
		
		int width = this.getWidth();
		x = (int)((0.5)*(this.getWidth() - dimension));
		//y = (int)((0.5)*(this.getHeight() - dimension));
		//System.out.println(this.getWidth());
		int height = this.getHeight();
		gr2.setPaint(Color.BLACK);
		Rectangle plot = new Rectangle(x, y, dimension, dimension);
		//plot.setFrame(x, y, dimension, dimension);
		
		//Rectangle plot = new Rectangle(100, 100, 500, 500);
		gr2.setPaint(Color.WHITE);
		gr2.fill(plot);
		
		gr2.setPaint(Color.BLUE);
		
		for(int i = x; i < x+dimension; i = i+OFFSET) {
			for(int j = y; j < y+dimension; j = j+OFFSET) {
				//float probability = random.nextFloat();
				float probability = probs[(i-x)/OFFSET][(j-y)/OFFSET];
				probability = 1/(float)Math.abs(Math.log(probability));
				//System.out.println("PROBABILITY: " + probability);
				gr2.setPaint(cg.getColor(probability));
				gr2.fill(new Ellipse2D.Double(i, j, OFFSET, OFFSET));
			}
		}
		
		count++;
		
	}
	
	public void setMatrix(float[][] probMatrix) {
		probs = probMatrix;
	}
	
	public void changeDimension(int x) {
		dimension = x;
	}
	
	public void clear() {
		
	}
	
	public void updateGUI(int i, int j, float probability) { // Update dot-plot area dynamically.
		//float probability = random.nextFloat();
		//plotArea.setBackground(cg.getColor(probability));
		
		
		
	}

	
	
	
	
	
}
