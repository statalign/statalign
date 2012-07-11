package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.util.ArrayList;

import javax.swing.JPanel;

import statalign.postprocess.plugins.LogLikelihoodTrace;
import statalign.postprocess.utils.LogLikelihoodTraceContainer;

/**
 * This class implements the graphical interface for showing the log-likelihood
 * plot.
 *
 * @author miklos, novak
 *
 */
public class LogLikelihoodTraceGUI extends JPanel {

	private static final Font LLT_FONT = new Font("Dialog", Font.PLAIN, 10);
	private static final long serialVersionUID = 1L;

	// int cornerx, cornery;
//	private JPanel parent;
	private LogLikelihoodTrace owner;

	/**
	 * Constructor to initialise the GUI for loglikelihood trace
	 *
	 * @param parent
	 *            The main panel
	 * @param owner
	 *            The logLikelihoodTrace postprocess handler
	 */
	public LogLikelihoodTraceGUI(JPanel parent, LogLikelihoodTrace owner) {
//		super((int) (panel.getWidth() / 6.6), panel.getHeight() / 13);
//		this.parent = parent;
		this.owner = owner;
//		setFont(new Font("Monospaced", Font.PLAIN, 10));
//		setEditable(false);
	}

	/**
	 * It updates the graphics of the panel
	 */
	@Override
	public void paintComponent(Graphics gr) {
		super.paintComponent(gr);
		
		final int border = 10;
		
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
		ArrayList<LogLikelihoodTraceContainer> list = owner.list;
		// finding the maximum and minimum
		double maxLik = -1000000000.0, minLik = 0.0;
		for (int i = 0; i < list.size(); i++) {
			double x = (list.get(i)).loglikelihood;
			if (x < minLik) {
				minLik = x;
			}
			if (x > maxLik) {
				maxLik = x;
			}
		}
		if (minLik > -0.1) {
			maxLik = 0.0;
		}
		gr.setFont(LLT_FONT);
		gr.drawString("" + ((int) maxLik), minX, 15+minY);
		gr.drawString("" + ((int) minLik), minX, maxHeight);
		// drawing the loglikelihood trace
		if (list.size() == 0) {
			return;
		}
		gr.setColor(new Color(221, 87, 20));
		double actual;
		double next = (list.get(0)).loglikelihood;
		boolean burnin = true;
		for (int i = 0; i < list.size() - 1; i++) {
			actual = next;
			next = (list.get(i + 1)).loglikelihood;
			if (burnin) {
				burnin = (list.get(i + 1)).burnin;
				if (!burnin) {
					gr.setColor(new Color(46, 87, 221));
				}
			}
			gr.drawLine(minX+(maxWidth-border) * i / 300 + 50,
					minY+(int) ((maxLik - actual) * (maxHeight-border) / (maxLik - minLik + 1.0)),
					minX+(maxWidth-border) * (i + 1) / 300 + 50,
					minY+(int) ((maxLik - next) * (maxHeight-border) / (maxLik - minLik + 1.0)));
			// i++;
		}
	}
	
}
