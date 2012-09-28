package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.ArrayList;

import javax.swing.JPanel;


import statalign.postprocess.plugins.VisualDistance;

public class DistanceGUI extends JPanel {

	/**
	 * 
	 * 
	 * @author Preeti Arunapuram
	 */
	
	private static final Font LLT_FONT = new Font("Dialog", Font.PLAIN, 10);
	private static final long serialVersionUID = 1L;
	
	//private JPanel parent;
	
	public String title;
	private VisualDistance owner;
	
	int distance = 1;
	
	public static final int OFFSET_X = 50;
	public static final int TITLE_X = 100;
	public static final int TITLE_Y = 30;

	public DistanceGUI(String title, VisualDistance owner) {
		this.title = title;
		this.owner = owner;
	}
	
	public void paintComponent(Graphics gr) {
		super.paintComponent(gr);
		
		Graphics2D g2 = (Graphics2D)gr;
		
		final int border = 10;
		
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, getWidth(), getHeight());
		g2.setColor(Color.black);
		
		int maxWidth = getWidth()-50-border;
		int maxHeight = getHeight()-2*border;
		int minX = border;
		int minY = border+30;
		g2.drawLine(50+minX, minY, 50+minX, maxHeight);
		g2.drawLine(50+minX, minY, 40+minX, 10+minY);
		g2.drawLine(50+minX, minY, 60+minX, 10+minY);
		g2.drawLine(50+minX, maxHeight, maxWidth + 50, maxHeight);
		g2.drawLine(maxWidth + 50, maxHeight, maxWidth + 40, maxHeight - 10);
		g2.drawLine(maxWidth + 50, maxHeight, maxWidth + 40, maxHeight + 10);
		ArrayList<Double> list = owner.distances;
		// finding the maximum and minimum
		double maxLik = 1.0, minLik = 0.0;
		for (int i = 0; i < list.size(); i++) {
			double x = list.get(i);
			if (x < minLik) {
				minLik = x;
			}
			if (x > maxLik) {
				maxLik = x;
			}
		}
		/*if (minLik > -0.1) {
			maxLik = 0.0;
		}*/
		g2.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		g2.rotate(Math.PI/2);
		//g2.rotate(-Math.PI*3/2);
		g2.setFont(new Font("SANS_SERIF", Font.PLAIN, 14));
		String yaxis = "Similarity (1st alignment, current)";
		g2.drawString(yaxis, (this.getHeight()-g2.getFontMetrics().stringWidth(yaxis))/2+5, -31);
		g2.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		//g2.rotate(Math.PI*3/2);
		g2.rotate(Math.PI*1.5);  
		g2.drawString("Sample", 700, 15+maxHeight);
		
		
		g2.drawString("" + ((int) maxLik), minX, 15+minY);
		g2.drawString("" + ((int) minLik), minX, maxHeight);
		// drawing the loglikelihood trace
		if (list.size() == 0) {
			g2.drawString("Waiting for data..", 100, 30);
			return;
		}
		
		g2.setColor(Color.BLACK);
		g2.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		g2.drawString(title, TITLE_X, TITLE_Y);
		g2.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		
		g2.setColor(Color.BLUE);
		double actual;
		double next = 1;
		for (int i = 0; i < list.size() - 1; i++) {
			actual = next;
			next = list.get(i);
			
			g2.drawLine(minX+(1000-border) * i * 2 / 300 + 50,
					minY+(int) ((maxLik - actual) * (maxHeight-border) / (maxLik - minLik + 1.0)),
					(minX+(1000-border) * (i + 1) * 2 / 300 + 50),
					minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)));
			// i++;
		}
	}
	
	/**
	 * This function tells the minimum size of the panel
	 */
	
	@Override
	public Dimension getMinimumSize(){
		return getPreferredSize();
	}

	/**
	 * This function tells the preferred size of the panel
	 */
	
	@Override
	public Dimension getPreferredSize() {
		//int maxWidth = this.getWidth() - 50 - border;
		return new Dimension(2*OFFSET_X + 7*(owner.mcmc.mcmcpars.cycles/owner.mcmc.mcmcpars.sampRate), TITLE_Y + 30);
//		if (alignment != null && alignment[0] != null) {
//			return new Dimension((alignment[0].length() + 3) * COLUMN_WIDTH + 6 * OFFSET_X, 100);
//		} else {
//			return new Dimension(0,100);
//		}
	}
	
}
	

