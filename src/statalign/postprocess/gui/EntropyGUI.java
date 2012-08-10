package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.util.ArrayList;

import javax.swing.JLabel;
import javax.swing.JPanel;

import statalign.postprocess.plugins.Entropy;
import statalign.postprocess.utils.EntropyContainer;

/**
 * This class implements the graphical interface for showing the observed, expected,
 * and individual sample entropies of a given alignment space.
 *
 * @author Preeti Arunapuram
 *
 */
public class EntropyGUI extends JPanel {

	private static final Font LLT_FONT = new Font("Dialog", Font.PLAIN, 10);
	private static final long serialVersionUID = 1L;

	// int cornerx, cornery;
//	private JPanel parent;
	
	public String title;
	private Entropy owner;

	private JLabel oe = new JLabel("Observed Entropy");
	private JLabel se = new JLabel("Sample Entropy");
	
	public static final int OFFSET_X = 50;
	
	public static final int TITLE_X = 100;
	public static final int TITLE_Y = 30;
	
	private int border = 10;
	
	/**
	 * Constructor to initialise the GUI for entropy
	 *
	 * @param parent
	 *            The main panel
	 * @param owner
	 *            The Entropy postprocess handler
	 */
	public EntropyGUI(String title, Entropy owner) {
//		super((int) (panel.getWidth() / 6.6), panel.getHeight() / 13);
//		this.parent = parent;
		this.title = title;
		this.owner = owner;
		//oe.setBackground(Color.BLUE);
		//se.setBackground(Color.RED);
		
		oe.setForeground(Color.BLUE);
		se.setForeground(Color.RED);
//		setFont(new Font("Monospaced", Font.PLAIN, 10));
//		setEditable(false);
	}

	/**
	 * It updates the graphics of the panel
	 */
	@Override
	public void paintComponent(Graphics gr) {
		//this.add(oe);
		//this.add(se);
		super.paintComponent(gr);
		
		border = 10;
		
		gr.setColor(Color.white);
		gr.fillRect(0, 0, getWidth(), getHeight());
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
		ArrayList<EntropyContainer> list = owner.entropyList;
		
		
		// finding the maximum and minimum
		double maxLik = 50.0, minLik = 0.0;
		for (int i = 0; i < list.size(); i++) {
			double x = (list.get(i)).obsEntropy;
			if (x < minLik) {
				minLik = x;
			}
			if (x + 20 > maxLik) {
				maxLik = x + 20;
			}
		}
		/*if (minLik > -0.1) {
			maxLik = 0.0;
		}*/
		gr.setFont(LLT_FONT);
		gr.drawString("" + ((int) maxLik), minX, 15+minY);
		gr.drawString("" + ((int) minLik), minX, maxHeight);
		// drawing the loglikelihood trace
		if (list.size() <= 1) {
			System.out.println("CYCLES AND SAMPLING RATE: " + owner.mcmc.mcmcpars.cycles + "; " + owner.mcmc.mcmcpars.sampRate);
			gr.drawString("Waiting for data..", TITLE_X, TITLE_Y);
			return;
		}
		
		System.out.println("CYCLES AND SAMPLING RATE HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH: " + owner.mcmc.mcmcpars.cycles + "; " + owner.mcmc.mcmcpars.sampRate);
		gr.setColor(Color.BLACK);
		gr.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		gr.drawString(title, TITLE_X, TITLE_Y);
		gr.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		
		this.add(oe);
		this.add(se);
		gr.setColor(Color.BLUE);
		double actual;
		double next = (list.get(0)).obsEntropy;
		for (int i = 0; i < list.size() - 1; i++) {
			actual = next;
			next = (list.get(i + 1)).obsEntropy;
			
			gr.drawLine(minX+(1000-border) * i * 2 / 300 + 50,
					minY+(int) ((maxLik - actual) * (maxHeight-border) / (maxLik - minLik + 1.0)),
					minX+(1000-border) * (i + 1) * 2 / 300 + 50,
					minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)));
			
			// i++;
			
			//gr.setColor(Color.BLUE);
			//gr.clearRect(minX+(maxWidth-border) * (i + 1) * 2 / 300 + 50, minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)), 100, 100);
			
			//gr.setColor(Color.WHITE);
			/*gr.drawString("Observed Entropy", (minX+(maxWidth-border) * (i + 1) * 2 / 300 + 50) + 10, 
					(minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0))) + 10);*/
			
			//JLabel oe = new JLabel("Observed Entropy");
			//oe.setForeground(Color.BLUE);
			//this.add(oe);
			oe.setLocation((minX+(1000-border) * (i + 1) * 2 / 300 + 50) + 10, 
					minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)));
		}
		
		//System.out.println("DRAWING LINES");
		/*gr.setColor(Color.GREEN);
		for (int i = 0; i < list.size() - 1; i++) {
			actual = next;
			next = (list.get(i + 1)).expEntropy;
			//System.out.println("DRAWING LINES");
			gr.drawLine(minX+(maxWidth-border) * i * 2 / 300 + 50,
					minY+(int) ((maxLik - actual) * (maxHeight-border) / (maxLik - minLik + 1.0)),
					minX+(maxWidth-border) * (i + 1) * 2 / 300 + 50,
					minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)));
			// i++;
		}*/
		
		gr.setColor(Color.RED);
		for (int i = 0; i < list.size() - 1; i++) {
			actual = next;
			next = (list.get(i + 1)).sampleEntropy;
			//System.out.println("DRAWING LINES");
			gr.drawLine(minX+(1000-border) * i * 2 / 300 + 50,
					minY+(int) ((maxLik - actual) * (maxHeight-border) / (maxLik - minLik + 1.0)),
					minX+(1000-border) * (i + 1) * 2 / 300 + 50,
					minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)));
			// i++;
			
			//JLabel se = new JLabel("Sample Entropy");
			//se.setForeground(Color.RED);
			//this.add(se);
			
			se.setLocation((minX+(1000-border) * (i + 1) * 2 / 300 + 50) + 10, 
					minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)));
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
