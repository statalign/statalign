package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
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
	
	int distance = -1;
	
	public static final int TITLE_X = 100;
	public static final int TITLE_Y = 30;

	public DistanceGUI(String title, VisualDistance owner) {
		this.title = title;
		this.owner = owner;
	}
	
	public void paintComponent(Graphics gr) {
		super.paintComponent(gr);
		
		final int border = 10;
		
		gr.setColor(Color.WHITE);
		gr.fillRect(19, 0, getWidth(), getHeight());
		gr.setColor(Color.black);
		
		int maxWidth = getWidth()-50-border;
		int maxHeight = getHeight()-2*border;
		int minX = border;
		int minY = border;
		gr.drawLine(50+minX, minY, 50+minX, maxHeight);
		gr.drawLine(50+minX, minY, 40+minX, 10+minY);
		gr.drawLine(50+minX, minY, 60+minX, 10+minY);
		gr.drawLine(50+minX, maxHeight, maxWidth + 50, maxHeight);
		gr.drawLine(maxWidth + 50, maxHeight, maxWidth + 40, maxHeight - 10);
		gr.drawLine(maxWidth + 50, maxHeight, maxWidth + 40, maxHeight + 10);
		ArrayList<Double> list = owner.distances;
		// finding the maximum and minimum
		double maxLik = 1.0, minLik = 0.0;
		/*for (int i = 0; i < list.size(); i++) {
			double x = (list.get(i)).obsEntropy;
			if (x < minLik) {
				minLik = x;
			}
			if (x > maxLik) {
				maxLik = x;
			}
		}
		if (minLik > -0.1) {
			maxLik = 0.0;
		}*/
		gr.setFont(LLT_FONT);
		gr.drawString("" + ((int) maxLik), minX, 15+minY);
		gr.drawString("" + ((int) minLik), minX, maxHeight);
		// drawing the loglikelihood trace
		if (list.size() == 0) {
			gr.drawString("Waiting for data..", 100, 30);
			return;
		}
		
		gr.setColor(Color.BLACK);
		gr.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		gr.drawString(title, TITLE_X, TITLE_Y);
		gr.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		
		gr.setColor(Color.BLUE);
		double actual;
		double next = 0;
		for (int i = 0; i < list.size() - 1; i++) {
			actual = next;
			next = list.get(i);
			
			gr.drawLine(minX+(maxWidth-border) * i * 2 / 300 + 50,
					minY+(int) ((maxLik - actual) * (maxHeight-border) / (maxLik - minLik + 1.0)),
					(minX+(maxWidth-border) * (i + 1) * 2 / 300 + 50),
					minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)));
			// i++;
		}
	}
	}
	

