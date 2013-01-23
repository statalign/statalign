package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.util.List;

import javax.swing.JPanel;

import statalign.postprocess.plugins.structalign.StructTrace;
import statalign.postprocess.utils.StructAlignTraceParameters;
import statalign.postprocess.utils.StructAlignTraceParameterGetters;
import statalign.postprocess.utils.StructAlignTraceParameterGetters.*;

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
	 * Constructor to initialise the GUI for StructAlign trace
	 *
	 * @param parent
	 *            The main panel
	 * @param owner
	 *            The StructTrace postprocess handler
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
		
		int maxX = getWidth()-50-border;
		int maxY = getHeight()-2*border;
		int minX = border;
		int minY = border;
		
		List<StructAlignTraceParameters> parameterHistory = owner.getParameterHistory();
		StructAlignTraceParameterGetters g = new StructAlignTraceParameterGetters(); 
		
		int plotSep = 30;
		if (parameterHistory.get(0).globalSigma) {
			paintParameter(gr,border,plotSep,minX,minY,
					maxX / 2 - plotSep, maxY,
					g.new Sigma2Getter(0), parameterHistory);
		}
		else {	
			int nSubplots = ( 1 + parameterHistory.get(0).sigma2.length );
			int subplotHeight = (int) ((double) maxY / (double) nSubplots);
			int i=0;
			for (i=0; i<(nSubplots-1); i++) {
				paintParameter(gr,border,plotSep,minX,minY+(i-1)*subplotHeight,
						maxX / 2 - plotSep, i*subplotHeight,
						g.new Sigma2Getter(i), parameterHistory);				
			}
			paintParameter(gr,border,plotSep,minX,minY+(i-1)*subplotHeight,
					maxX / 2 - plotSep, i*subplotHeight,
					g.new Sigma2HGetter(), parameterHistory);
			
		}
		int nSubplots = 3; // nu, tau, epsilon
		int subplotHeight = (int) ((double) maxY / (double) nSubplots);
		paintParameter(gr,border,plotSep,minX+maxX/2 + plotSep,minY+0*subplotHeight,
			maxX - border, 1*subplotHeight,
			g.new NuGetter(), parameterHistory);				
		paintParameter(gr,border,plotSep,minX+maxX/2 + plotSep,minY+1*subplotHeight,
			maxX - border, 2*subplotHeight,
			g.new TauGetter(), parameterHistory);
		paintParameter(gr,border,plotSep,minX+maxX/2 + plotSep,minY+2*subplotHeight,
			maxX - border, 3*subplotHeight,
			g.new EpsilonGetter(), parameterHistory);
	
		
	}
	
	private void paintParameter(Graphics gr, final int border, int plotSep,
			int minX, int minY, int maxX, int maxY,
			ParameterGetter getter, List<StructAlignTraceParameters> parameterHistory) {
		
		// y-axis
		gr.drawLine(20+minX,minY,20+minX,maxY);
		// x-axis
		gr.drawLine(20+minX,minY,50+maxX,minY);
		// right border
		gr.drawLine(50+maxX,minY,50+maxX,maxY);
		// top border
		gr.drawLine(20+minX,maxY,50+maxX,maxY);
		
		//gr.drawLine(20+minX, minY, 10+minX, 10+minY);
		//gr.drawLine(20+minX, minY, 30+minX, 10+minY);
		
		//gr.drawLine(maxX + 50, maxY, maxX + 40, maxY - 10);
		//gr.drawLine(maxX + 50, maxY, maxX + 40, maxY + 10);
		// finding the maximum and minimum
		double maxParam = 0.0, minParam = 10000000.0;
		for (int i = 0; i < parameterHistory.size(); i++) {
			double param = getter.getParameter(parameterHistory.get(i));
			if (param < minParam) {
				minParam = param;
			}
			if (param > maxParam) {
				maxParam = param;
			}
		}
		gr.setFont(STR_FONT);
		gr.drawString("" + ((int) maxParam), minX, 15+minY);
		gr.drawString("" + ((int) minParam), minX, 10 + maxY);
		// TODO replace the (int) casts with sprintf to 
		// 3sf.
		
		// drawing the trace
		if (parameterHistory.size() == 0) {
			return;
		}
		gr.setColor(new Color(221, 87, 20));
		double current;
		double next = getter.getParameter(parameterHistory.get(0));
		boolean burnin = true;
		for (int i = 0; i < parameterHistory.size() - 1; i++) {
			current = next;
			next = getter.getParameter(parameterHistory.get(i + 1));
			if (burnin) {
				burnin = (parameterHistory.get(i + 1)).burnin;
				if (!burnin) {
					gr.setColor(new Color(46, 87, 221));
				}
			}
			gr.drawLine(minX+(maxX-minX) * i / 300 + 20,
					minY+(int) ((maxParam - current) * (maxY-minY) / (maxParam - minParam + 1.0)),
					minX+(maxX-minX) * (i + 1) / 300 + 20,
					minY+(int) ((maxParam - next) * (maxY-minY) / (maxParam - minParam + 1.0)));
		}
		gr.setColor(Color.black);
	}
}
