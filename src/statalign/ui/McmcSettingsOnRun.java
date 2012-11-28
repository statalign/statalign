package statalign.ui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import statalign.base.AutomateParameters;
import statalign.base.MCMCPars;
import statalign.base.AutomateParameters;

/**
 * 
 * This is the dialog window where users can set MCMC parameters.
 * 
 * @author miklos, novak
 *
 */
public class McmcSettingsOnRun extends JDialog implements ActionListener, KeyListener {
	private static final long serialVersionUID = 1L;

	private MCMCPars pars;
//	boolean toRun = false;
	
	JTextField burnIn = new JTextField(10);
	JTextField cycles = new JTextField(10);
	JTextField sampRate = new JTextField(10);
	JTextField seed = new JTextField(10);
	JCheckBox automateStepRate = new JCheckBox("Automate (Slow)",true);
	JCheckBox automateNumberOfSamples = new JCheckBox("Automate (RNA only)",false);
	JCheckBox automateBurnIn = new JCheckBox("Automate",true);
//	private JTextField outFile = new JTextField(15)
	private MainFrame owner;
	
	McmcSettingsOnRun(MainFrame owner) {
		super(owner, "MCMC Parameters", true);
		this.owner = owner;
		pars = owner.manager.inputData.pars;
		Container cp = getContentPane();
//		cp.setLayout(new BoxLayout(cp,BoxLayout.Y_AXIS));
		cp.setLayout(new BorderLayout());
		Box bigBox = Box.createVerticalBox();
		JPanel pan = new JPanel();
		GridLayout l = new GridLayout(4,4);
		l.setHgap(5);
		l.setVgap(5);
		pan.setLayout(l);
		pan.add(new JLabel("Burn-in cycles:"));
		burnIn.addKeyListener(this);
		pan.add(burnIn);
		burnIn.setEnabled(false);
		
		//pan.add(new JLabel("Autometic"));
		automateBurnIn.setActionCommand("burnin");
		automateBurnIn.addKeyListener(this);
		automateBurnIn.addActionListener(this);
		pan.add(automateBurnIn);
		
		pan.add(new JLabel("Cycles after burn-in:"));
		cycles.addKeyListener(this);
		pan.add(cycles);
		cycles.setEnabled(true);
		
		//pan.add(new JLabel("Autometic"));
		automateNumberOfSamples.setActionCommand("numsam");
		automateNumberOfSamples.addKeyListener(this);
		automateNumberOfSamples.addActionListener(this);
		pan.add(automateNumberOfSamples);
		
		pan.add(new JLabel("Sampling rate:"));
		sampRate.addKeyListener(this);
		pan.add(sampRate);
		sampRate.setEnabled(false);
		
		
		//pan.add(new JLabel("Autometic"));
		automateStepRate.setActionCommand("steprate");
		automateStepRate.addKeyListener(this);
		automateStepRate.addActionListener(this);
		pan.add(automateStepRate);
		
		pan.add(new JLabel("Seed:"));
		seed.addKeyListener(this);
		pan.add(seed);
		
		
		
//		pan.add(new JLabel("Output file:"));
//		pan.add(outFile);
		bigBox.add(pan);
		Box box = Box.createHorizontalBox();
		JButton butt;
		box.add(butt=new JButton("Run!"));
		butt.addActionListener(this);
		getRootPane().setDefaultButton(butt);
		box.add(Box.createHorizontalStrut(100));
		box.add(butt=new JButton("Cancel"));
		butt.addActionListener(this);
		bigBox.add(box);
		cp.add(bigBox, SwingConstants.CENTER);
		cp.add(Box.createHorizontalStrut(20), BorderLayout.LINE_START);
		cp.add(Box.createHorizontalStrut(20), BorderLayout.LINE_END);
		cp.add(Box.createVerticalStrut(15), BorderLayout.PAGE_START);
		cp.add(Box.createVerticalStrut(15), BorderLayout.PAGE_END);
		addKeyListener(this);
		pack();
//		bigBox.setMaximumSize(bigBox.getSize());
//		setSize(getWidth()+30,getHeight()+30);
	}
	
	void display(Component c) {
		burnIn.setText(Integer.toString(pars.burnIn));
		cycles.setText(Integer.toString(pars.cycles));
		sampRate.setText(Integer.toString(pars.sampRate));
		seed.setText(Integer.toString((int) pars.seed));
//		outFile.setText(sp.outFile);
		setLocationRelativeTo(c);
//		pack();
		setVisible(true);
	}

	/**
	 * 
	 * This is inherited from the ActionListener interface.
	 * When we close the dialog, it updates the MCMC parameters.
	 * 
	 */
	public void actionPerformed(ActionEvent ev) {
		if(ev.getActionCommand() == "numsam"){
			if(automateNumberOfSamples.isSelected()){
				cycles.setEnabled(false);
				
				owner.mcmcSettingsDlg.automateNumberOfSamples.setSelected(true);
				owner.mcmcSettingsDlg.cycles.setEnabled(false);
			}
			else{
				cycles.setEnabled(true);
				
				owner.mcmcSettingsDlg.automateNumberOfSamples.setSelected(false);
				owner.mcmcSettingsDlg.cycles.setEnabled(true);
			}
		}
		if(ev.getActionCommand() == "steprate"){
			if(automateStepRate.isSelected()){
				sampRate.setEnabled(false);
				
				owner.mcmcSettingsDlg.automateStepRate.setSelected(true);
				owner.mcmcSettingsDlg.sampRate.setEnabled(false);
			}
			else{
				sampRate.setEnabled(true);
				
				owner.mcmcSettingsDlg.automateStepRate.setSelected(false);
				owner.mcmcSettingsDlg.sampRate.setEnabled(true);
			}
		}
		if(ev.getActionCommand() == "burnin"){
			if(automateBurnIn.isSelected()){
				burnIn.setEnabled(false);
				
				owner.mcmcSettingsDlg.automateBurnIn.setSelected(true);
				owner.mcmcSettingsDlg.burnIn.setEnabled(false);
			}
			else{
				burnIn.setEnabled(true);
				
				owner.mcmcSettingsDlg.automateBurnIn.setSelected(false);
				owner.mcmcSettingsDlg.burnIn.setEnabled(true);
			}
		}
		
		try{
			if(ev.getActionCommand() == "Run!") {
				pars.burnIn = Integer.parseInt(burnIn.getText());
				pars.cycles = Integer.parseInt(cycles.getText());
				pars.sampRate = Integer.parseInt(sampRate.getText());
				pars.seed = Integer.parseInt(seed.getText());
				AutomateParameters.setAutomateStepRate(automateStepRate.isSelected());
				AutomateParameters.setAutomateNumberOfSamples(automateNumberOfSamples.isSelected());
				AutomateParameters.setAutomateBurnIn(automateBurnIn.isSelected());
//				sp.outFile = outFile.getText();
//				toRun = true;
				setVisible(false);
				
				owner.disableAllButtons();
				owner.start();
			}
			if(ev.getActionCommand() == "Cancel") {
				setVisible(false);
			}
//				toRun = false;
			
			
		}
		catch(NumberFormatException e){
			new ErrorMessage(owner,"Wrong format, "+e.getLocalizedMessage(),false);
		}
	}

	public void keyPressed(KeyEvent e) {}
	public void keyTyped(KeyEvent e) {}

	public void keyReleased(KeyEvent e) {
		if(e.getKeyCode() == KeyEvent.VK_ESCAPE) {
			setVisible(false);
		}
	}
	
}
