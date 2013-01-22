package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.util.List;

import javax.swing.JPanel;

import statalign.postprocess.plugins.structalign.StructTrace;
import statalign.postprocess.utils.StructAlignTraceParameters;

/**
 * This class implements the graphical interface for showing the log-likelihood
 * plot.
 *
 * @author miklos, novak
 *
 */
public class StructAlignTraceGUI extends JPanel {

	private static final Font STR_FONT = new Font("Dialog", Font.PLAIN, 10);
	private static final long serialVersionUID = 1L;

	// int cornerx, cornery;
//	private JPanel parent;
	private StructTrace owner;

	/**
	 * Constructor to initialise the GUI for loglikelihood trace
	 *
	 * @param parent
	 *            The main panel
	 * @param owner
	 *            The logLikelihoodTrace postprocess handler
	 */
	public StructAlignTraceGUI(JPanel parent, StructTrace owner) {
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
		List<StructAlignTraceParameters> parameterHistory = owner.getParameterHistory();
		// finding the maximum and minimum
		double maxTau = 0.0, minTau = 10000000.0;
		for (int i = 0; i < parameterHistory.size(); i++) {
			double tau = (parameterHistory.get(i)).tau;
			if (tau < minTau) {
				minTau = tau;
			}
			if (tau > maxTau) {
				maxTau = tau;
			}
		}
		gr.setFont(STR_FONT);
		gr.drawString("" + ((int) maxTau), minX, 15+minY);
		gr.drawString("" + ((int) minTau), minX, maxHeight);
		// drawing the trace
		if (parameterHistory.size() == 0) {
			return;
		}
		gr.setColor(new Color(221, 87, 20));
		double actual;
		double next = (parameterHistory.get(0)).tau;
		boolean burnin = true;
		for (int i = 0; i < parameterHistory.size() - 1; i++) {
			actual = next;
			next = (parameterHistory.get(i + 1)).tau;
			if (burnin) {
				burnin = (parameterHistory.get(i + 1)).burnin;
				if (!burnin) {
					gr.setColor(new Color(46, 87, 221));
				}
			}
			gr.drawLine(minX+(maxWidth-border) * i / 300 + 50,
					minY+(int) ((maxTau - actual) * (maxHeight-border) / (maxTau - minTau + 1.0)),
					minX+(maxWidth-border) * (i + 1) / 300 + 50,
					minY+(int) ((maxTau - next) * (maxHeight-border) / (maxTau - minTau + 1.0)));
			// i++;
		}
	}
	
}
