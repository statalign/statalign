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

import statalign.base.AutomateParamSettings;
import statalign.base.MCMCPars;

/**
 * 
 * This is the dialog window where users can set MCMC parameters.
 * 
 * @author miklos, novak
 *
 */
public class McmcSettingsDlg extends JDialog implements ActionListener, KeyListener {
	private static final long serialVersionUID = 1L;

	private MCMCPars pars;
//	boolean toRun = false;
	
	JTextField burnIn = new JTextField(10);
	JTextField cycles = new JTextField(10);
	JTextField sampRate = new JTextField(10);
	JTextField seed = new JTextField(10);
	JCheckBox automateStepRate = new JCheckBox("Automate (Slow)");
	JCheckBox automateNumberOfSamples = new JCheckBox("Automate (RNA only)");
	JCheckBox automateBurnIn = new JCheckBox("Automate");
//	private JTextField outFile = new JTextField(15)
	private MainFrame owner;
	
	McmcSettingsDlg(MainFrame owner) {
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
		
		//pan.add(new JLabel("Automatic"));
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
		bigBox.add(Box.createVerticalStrut(20));
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
		AutomateParamSettings autoPars = pars.autoParamSettings;
		automateBurnIn.setSelected(autoPars.automateBurnIn);
		automateNumberOfSamples.setSelected(autoPars.automateNumberOfSamplesToTake);
		automateStepRate.setSelected(autoPars.automateSamplingRate);
		
		updateEnabled();
		setLocationRelativeTo(c);
		setVisible(true);
	}

	/**
	 * 
	 * This is inherited from the ActionListener interface.
	 * When we close the dialog, it updates the MCMC parameters.
	 * 
	 */
	@Override
	public void actionPerformed(ActionEvent ev) {
		if(ev.getActionCommand() == "numsam" || ev.getActionCommand() == "steprate" ||
				ev.getActionCommand() == "burnin") {
					
			updateEnabled();
			
		} else if(ev.getActionCommand() == "Run!") {
			try {
				pars.burnIn = Integer.parseInt(burnIn.getText());
				pars.cycles = Integer.parseInt(cycles.getText());
				pars.sampRate = Integer.parseInt(sampRate.getText());
				pars.seed = Integer.parseInt(seed.getText());
				AutomateParamSettings autoPar = pars.autoParamSettings;
				autoPar.automateSamplingRate = automateStepRate.isSelected();
				autoPar.automateNumberOfSamplesToTake = automateNumberOfSamples.isSelected();
				autoPar.automateBurnIn = automateBurnIn.isSelected();
//				sp.outFile = outFile.getText();
//				toRun = true;
				setVisible(false);
				
				owner.disableAllButtons();
				owner.start();
			} catch(NumberFormatException e){
				ErrorMessage.showPane(owner, "Wrong format, "+e.getLocalizedMessage(), false);
			}
			
		} else if(ev.getActionCommand() == "Cancel") {
			setVisible(false);
		}
			
	}

	private void updateEnabled() {
		cycles.setEnabled(!automateNumberOfSamples.isSelected());
		sampRate.setEnabled(!automateStepRate.isSelected());
		burnIn.setEnabled(!automateBurnIn.isSelected());
	}

	@Override
	public void keyPressed(KeyEvent e) {}
	@Override
	public void keyTyped(KeyEvent e) {}

	@Override
	public void keyReleased(KeyEvent e) {
		if(e.getKeyCode() == KeyEvent.VK_ESCAPE) {
			setVisible(false);
		}
	}
	
}
