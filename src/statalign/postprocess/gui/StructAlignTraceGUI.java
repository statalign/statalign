package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.BasicStroke;
import java.awt.geom.Line2D;
import java.util.List;

import javax.swing.JPanel;

import statalign.postprocess.plugins.structalign.StructAlignTraceParameters;
import statalign.postprocess.plugins.structalign.StructTrace;
import statalign.postprocess.plugins.structalign.StructAlignTraceParameters.*;

/**
 * This class implements the graphical interface for showing the trace
 * plots for the StructAlign parameters.
 *
 * @author herman
 *
 */
public class StructAlignTraceGUI extends JPanel {

	private static final Font STR_FONT = new Font("Dialog", Font.PLAIN, 10);
	private static final long serialVersionUID = 1L;


	private StructTrace owner;
	private int burninLength;

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
		this.burninLength = owner.burninLength;
//		setFont(new Font("Monospaced", Font.PLAIN, 10));
//		setEditable(false);
	}

	/**
	 * It updates the graphics of the panel
	 */
	@Override
	public void paintComponent(Graphics gra) {
		super.paintComponent(gra);
		
		final int border = 10;
		
		gra.setColor(Color.white);
		gra.fillRect(0, 0, getWidth(), getHeight());
		gra.setColor(Color.black);
		
		Graphics2D gr = (Graphics2D) gra;
		
		int maxX = getWidth()-50-border;
		int maxY = getHeight()-2*border;
		int minX = border;
		int minY = border;
		
		List<StructAlignTraceParameters> parameterHistory = owner.getParameterHistory();
			
		int plotSep = 30;
		List<PlottableParameter> plottableParameters = parameterHistory.get(0).plottableParameters;
		int nPlottableParams = plottableParameters.size();
		int nLeft = 0, nRight = 0;
		for (int i=0; i<nPlottableParams; i++) {
			int plotSide = plottableParameters.get(i).plotSide;
			if (plotSide == 0) {
				nLeft++;
			}
			else {
				nRight++;
			}
		}
		int iLeft = 0, iRight = 0;
		for (int i=0; i<nPlottableParams; i++) {
			int plotSide = plottableParameters.get(i).plotSide;
			if (plotSide == 0) {
				int nSubplots = nLeft;
				double subplotHeight = (double) maxY / (double) nSubplots;
				paintParameter(gr,border,plotSep,minX,minY+iLeft*subplotHeight,
						maxX / 2 - plotSep, (iLeft+1)*subplotHeight,
						i, parameterHistory);
				iLeft++;
			}
			else {
				int nSubplots = nRight;
				double subplotHeight = (double) maxY / (double) nSubplots;
				paintParameter(gr,border,plotSep,minX+maxX/2 + plotSep,minY+iRight*subplotHeight,
						maxX - border, (iRight+1)*subplotHeight,
						i, parameterHistory);
				iRight++;
			}
		}	
	}
	
	private void paintParameter(Graphics2D gr, final int border, double plotSep,
			double minX, double minY, double maxX, double maxY,
			int parameterIndex, List<StructAlignTraceParameters> parameterHistory) {
		
		Shape line;
		
		double textShift = 30;
		// y-axis
		line = new Line2D.Double(textShift+minX,minY,textShift+minX,maxY);
		gr.draw(line);
		// x-axis
		line = new Line2D.Double(textShift+minX,minY,50+maxX,minY);
		gr.draw(line);
		// right border
		line = new Line2D.Double(50+maxX,minY,50+maxX,maxY);
		gr.draw(line);
		// top border
		line = new Line2D.Double(textShift+minX,maxY,50+maxX,maxY);
		gr.draw(line);
		
		//gr.drawLine(20+minX, minY, 10+minX, 10+minY);
		//gr.drawLine(20+minX, minY, 30+minX, 10+minY);
		
		//gr.drawLine(maxX + 50, maxY, maxX + 40, maxY - 10);
		//gr.drawLine(maxX + 50, maxY, maxX + 40, maxY + 10);
		// finding the maximum and minimum
		double maxParam = 0.0, minParam = 10000000.0;
		//for (int i = Math.min(owner.MAX_HISTORY_SIZE/4,parameterHistory.size() / 2); i < parameterHistory.size(); i++) {
		for (int i = 0; i < parameterHistory.size(); i++) {
			double param = (Double) parameterHistory.get(i).plottableParameters.get(parameterIndex).value;
			if (param < minParam) {
				minParam = param;
			}
			if (param > maxParam) {
				maxParam = param;
			}
		}
		gr.setFont(STR_FONT);
		Color burninColor = new Color(251, 87, 20);
		Color nonBurninColor = new Color(46, 87, 221);
		
		gr.drawString("" + String.format("%.1f",maxParam), (int) minX, 5 + (int)minY);
		gr.drawString("" + String.format("%.1f",minParam), (int) minX, 5 + (int)maxY);
		gr.setFont(new Font("Dialog", Font.BOLD, 12));
		String paramName = parameterHistory.get(0).plottableParameters.get(parameterIndex).name;
		if (parameterHistory.get(parameterHistory.size()-1).plottableParameters.get(parameterIndex).fixedToParent) {
			gr.setColor(burninColor);
		}
		gr.drawString(paramName, (int) minX + 3, (int)((minY+maxY)/2));
		gr.setFont(STR_FONT);
				
		// drawing the trace
		if (parameterHistory.size() == 0) {
			return;
		}
				
		gr.setColor(burninColor);
		//boolean finishedBurnin = (parameterHistory.get(parameterHistory.size()-1)).burnin;
		int startFrom = 0;
//		if (finishedBurnin) {
//			startFrom = parameterHistory.size() - 1 - burninLength; // need -1?
//			if (startFrom > burninLength) {
//				startFrom = burninLength; // i.e. start at first non-burnin step?
//			} 
//		}
		double current;
		double next = (Double) parameterHistory.get(startFrom).plottableParameters.get(parameterIndex).value;
//		if (next < minParam) {
//			next = minParam;
//		}
//		if (next > maxParam) {
//			next = maxParam;
//		}
		boolean burnin = true;
		for (int i = startFrom; i < parameterHistory.size() - 1; i++) {
			current = next;
			next = (Double) parameterHistory.get(i+1).plottableParameters.get(parameterIndex).value;
//			if (next < minParam) {
//				next = minParam;
//			}
//			if (next > maxParam) {
//				next = maxParam;
//			}
			if (burnin) {
				burnin = (parameterHistory.get(i + 1)).burnin;
			}
			boolean wasProposed = parameterHistory.get(i+1).plottableParameters.get(parameterIndex).wasProposed;
			if (!wasProposed) {
                gr.setStroke(new BasicStroke(0));
				gr.setColor(Color.black);
			}
			if (wasProposed) {
                gr.setStroke(new BasicStroke(2));
				if (burnin) {
					gr.setColor(burninColor);
				}
				else {
					gr.setColor(nonBurninColor);
				}
			}
				line = new Line2D.Double(minX+(maxX-minX) * i / parameterHistory.size() + textShift, 
            		minY+((maxParam - current) * (maxY-minY) / (maxParam - minParam)),
					minX+(maxX-minX) * (i + 1) / parameterHistory.size() + textShift,
					minY+ ((maxParam - next) * (maxY-minY) / (maxParam - minParam)));
				gr.draw(line);
			//}
//			gr.drawLine(minX+(maxX-minX) * i / 300 + 20,
//					minY+(int) ((maxParam - current) * (maxY-minY) / (maxParam - minParam + 1.0)),
//					minX+(maxX-minX) * (i + 1) / 300 + 20,
//					minY+(int) ((maxParam - next) * (maxY-minY) / (maxParam - minParam + 1.0)));
		}
		gr.setColor(Color.black);
		gr.setStroke(new BasicStroke());
		
		gr.setFont(new Font("Dialog", Font.BOLD, 10));
		double acceptanceRate = (Double) parameterHistory.get(parameterHistory.size()-1).plottableParameters.get(parameterIndex).acceptanceRate;
		gr.drawString("α = " + String.format("%.2f",acceptanceRate), 5 + (int)textShift + (int) minX, 15 + (int)minY);
		gr.setFont(STR_FONT);
	}
}
