package statalign.model.ext.plugins.structalign;

import statalign.model.ext.ModelExtension;
import statalign.model.ext.McmcMove;
import statalign.model.ext.plugins.StructAlign;

public abstract class StructAlignMcmcMove extends McmcMove {

	StructAlign owner;
	@Override
	public ModelExtension getOwner() {
		return owner;
	}
	
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
}
