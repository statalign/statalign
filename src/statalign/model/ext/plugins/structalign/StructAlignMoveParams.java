package statalign.model.ext.plugins.structalign;

public class StructAlignMoveParams {

	/** 
	 * <tt>true</tt> if a trace plot for this move will be displayed
	 * on the GUI panel.
	 */
	private boolean plottable = false;
	public boolean isPlottable() {
		return plottable;
	}
	public void setPlottable() {
		plottable = true;
	}
	public void unsetPlottable() {
		plottable = false;
	}
	/** 
	 * <tt>0</tt> if the trace plot should be plotted on the left side
	 * of the panel, and <tt>1</tt> if it should be plotted on the right. 
	 */
	private int plotSide;
	public int plotSide() {
		return plotSide;
	}
	public void setPlotSide(int s) {
		plotSide = s;
	}
	private boolean fixedToParent;
	public boolean fixedToParent() {
		return fixedToParent;
	}
	public void setFixedToParent(boolean b) {
		fixedToParent = b;
	}
}
