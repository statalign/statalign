package statalign.postprocess.plugins;

//import java.awt.*;

import java.awt.BorderLayout;
import java.util.ArrayList;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.base.InputData;
import statalign.base.State;
import statalign.postprocess.Postprocess;
import statalign.postprocess.gui.EntropyGUI;
import statalign.postprocess.utils.EntropyContainer;

/**
 * This is the implementation of the postprocess that handles entropy.
 *
 * @author Preeti Arunapuram
 *
 */
public class Entropy extends statalign.postprocess.Postprocess{

	public String title;
	
	JPanel pan = new JPanel(new BorderLayout());
	public ArrayList<EntropyContainer> entropyList;
	PPFold ppFold;
	
	public static boolean allowed = false;
	
	
	private EntropyGUI gui;

	/**
	 * It constructs a postprocess handler that is  screenable (=can be shown on the GUI),
	 * outputable (= can be written into the log file) and postprocessable
	 * (= can generate its own output file).
	 *
	 * The default setting is that it writes loglikelihoods both into the logfile and its own
	 * output file.
	 */
	public Entropy(){
		screenable = true;
		outputable = true;
//		postprocessable = true;
//		postprocessWrite = true;
		sampling = true;
		rnaAssociated = true;
		selected = false;
	}

	/**
	 * Returns the name of the postprocess, 'Loglikelihood'.
	 */
	@Override
	public String getTabName(){
		//return null;
		return "Entropy";
	}

	/**
	 * Returns the icon generated from the picture in file icons/loglikelihood1.gif.
	 */
	@Override
	public Icon getIcon(){
		return new ImageIcon(ClassLoader.getSystemResource("icons/loglikelihood1.gif"));
//		return new ImageIcon("icons/loglikelihood1.gif");
	}

	/**
	 * It generates a panel and returns it.
	 */
	@Override
	public JPanel getJPanel(){
		return pan;
	}
	
	@Override
	public void reloadPanel() {
		pan = new JPanel(new BorderLayout());
	}

	/**
	 * Returns the tip of the tabulated panel 'Log-likelihood trace'
	 * (shown when mouse is moved over the panel)
	 */
	@Override
	public String getTip(){
		return "Entropy";
	}

	@Override
	public double getTabOrder() {
		return 10.0d;
	}
	
	@Override
	public String getFileExtension() {
		return "entr";
	}
	
	@Override
	public String[] getDependencies() {
		return new String[] { "statalign.postprocess.plugins.PPFold" };
	}
	
	@Override
	public void refToDependencies(Postprocess[] plugins) {
		ppFold = (PPFold)plugins[0];
	}

	/**
	 * Writes the loglikelihoods into a logfile, if sampling mode is switched on.
	 */
	@Override
	public void newSample(State state, int no, int total) {
		if(sampling){
			if(ppFold == null) {return;}
			
			//System.out.println("Observed Entropy Object: " + ppFold.entropyObs);
			//System.out.println(entropyList);
			//entropyList.add(new EntropyContainer(ppFold.fuzzyResultObs.entropyVal, ppFold.fuzzyResultExp.entropyVal));
			
			if(allowed) {
				//System.out.println("Observed Entropy Object: " + ppFold.entropyObs);
				entropyList.add(new EntropyContainer(ppFold.entropyObs, ppFold.entropyExp, ppFold.entropySample));
				
			}
			
			if(show) {
				gui.repaint();
			}
			
		}
		else{
			//System.out.println("Not sampling loglikelihood!!!");
		}
	}

	/**
	 * Initialises the graphical interface.
	 */
	@Override
	public void beforeFirstSample(InputData input) {
		
		entropyList = new ArrayList<EntropyContainer>();
		
		if(show) {
			title = input.title;
			pan.removeAll();
			gui = new EntropyGUI(title, this);
			
			JScrollPane scroll = new JScrollPane();
			scroll.setViewportView(gui);
			
			pan.add(scroll);
			if(pan.getParent() != null)
			{
				pan.getParent().getParent().getParent().validate();
			}
		}
		
		
	}

	/**
	 * It writes the sampled loglikelihood values into a separate file.
	 */
	@Override
	public void afterLastSample() {
	}

	/**
	 * It switches on and off the sampling mode.
	 */
	@Override
	public void setSampling(boolean enabled) {
		sampling = enabled;

	}

}