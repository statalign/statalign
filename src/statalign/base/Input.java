package statalign.base;

//import java.awt.*;

import java.awt.BorderLayout;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;

import statalign.postprocess.gui.InputGUI;

/**
 * 
 * This class is an extension of the Postprocess class and is for showing the
 * actual input sequences.
 * 
 * @author miklos, novak
 *
 */
public class Input extends statalign.postprocess.Postprocess{

	JPanel pan;
	MainManager manager;
	
	/**
	 * this is the graphical interface
	 */
	public InputGUI inputgui;
	
	/**
	 * A trivial constructor. Only sets the MainManager of this object.
	 * @param manager the MainManager via that the object could communicate and can get
	 * information.
	 */
	public Input(MainManager manager){
		this.manager = manager;
		screenable = true;
	}
	
	@Override
	public JPanel getJPanel() {
		pan = new JPanel(new BorderLayout());
		inputgui = new InputGUI(manager); 
		pan.add(inputgui);
		return pan;
	}

	@Override
	public Icon getIcon() {
		return new ImageIcon(ClassLoader.getSystemResource("icons/sequences.gif"));
	}

	@Override
	public String getTabName() {
		return "Sequences";
	}

	@Override
	public String getTip() {
		return "Input sequences";
	}

	@Override
	public double getTabOrder() {
		return 1.0d;
	}

	/**
	 * This is an empty function.
	 * 
	 * Indeed we do not have to do anything with the input data when the MCMC makes a new step.
	 */
	public void newStep(){
		
	}

	/**
	 * This is an empty function.
	 * 
	 * Indeed we do not have to do anything with the input data before the first step in the MCMC run.
	 */
	public void beforeFirstSample() {
		//pan.add(new InputGUI(pan, this, mcmc.tree.rowSequences()));
		
	}

	/**
	 * This is an empty function.
	 * 
	 * Indeed we do not have to do anything with the input data after the last step in the MCMC run.
	 */
	@Override
	public void afterLastSample() {
		// Empty
	}

	/**
	 * This is an empty function.
	 * 
	 * Never used but has to be implemented as the class extends the Prostrocess class 
	 */
	@Override
	public void setSampling(boolean enabled) {
		// void function...
		
	}

}
