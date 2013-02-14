package statalign.ui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.File;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import statalign.postprocess.Postprocess;
import statalign.postprocess.utils.RNAalifold;

/**
 * 
 * This is the dialog window where users can set MCMC parameters.
 * 
 * @author miklos, novak
 *
 */
public class RNASettingsDlg extends JDialog implements ActionListener, ChangeListener, KeyListener {
	private static final long serialVersionUID = 1L;
	

	JPanel rnaalifoldOptions = new JPanel();
	JCheckBox useSamplingAndAveragingButton = new JCheckBox("Perform sampling and averaging prediction (PPfold).");
	JCheckBox useSamplingAndAveragingRNAalifoldButton = new JCheckBox("Perform sampling and averaging prediction (RNAalifold).");
	JCheckBox fuzzyNucleotidePredictionAndEntropy = new JCheckBox("Perform fuzzy alignment prediction (PPfold).");
	JCheckBox consensusEvolutionPrediction = new JCheckBox("Perform consensus alignment prediction (PPfold).");
	
	private JFileChooser fileChooser = new JFileChooser();
	private JTextField executableField = new JTextField("");
	private JButton executableButton = new JButton("Browse");
	
	private SpinnerNumberModel temperatureModel = new SpinnerNumberModel(37.0, -273.0, 500.0, 1.0);
	private JSpinner temperatureSpinner = new JSpinner(temperatureModel);
	
	JRadioButton linearButton = new JRadioButton("Linear");
	JRadioButton circularButton = new JRadioButton("Circular");
	//private JTextField cycles = new JTextField(10);	
	private SpinnerNumberModel covarianceModel = new SpinnerNumberModel(1.0, 0.0, 10000, 0.5);
	private JSpinner covarianceSpinner = new JSpinner(covarianceModel);
	private SpinnerNumberModel nonCompatibleModel = new SpinnerNumberModel(1.0, 0.0, 10000, 0.5);
	private JSpinner nonCompatibleSpinner = new JSpinner(nonCompatibleModel);
	private JFrame owner;
	
	private Preferences prefs;
	
	RNASettingsDlg(JFrame owner) {
		super(owner, "RNA options", true);
		this.owner = owner;
		prefs = Preferences.userNodeForPackage(this.getClass());
		
		//pars = owner.manager.inputData.pars;
		//Container cp = getContentPane();
		setLayout(new BorderLayout());
		Box bigBox = Box.createVerticalBox();
		GridBagLayout l = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 0.5;
		rnaalifoldOptions.setLayout(l);
		
		JPanel optionsPanel = new JPanel();
		GridLayout optionLayout = new GridLayout(3,1);
		optionsPanel.setLayout(optionLayout);
		//JCheckBox useSamplingAndAveragingButton = new JCheckBox("Use sampling and averaging prediction (PPfold).");
		//JCheckBox useSamplingAndAveragingRNAalifoldButton = new JCheckBox("Use sampling and averaging prediction (RNAalifold).");
		//JCheckBox fuzzyNucleotidePredictionAndEntropy
		useSamplingAndAveragingButton.setSelected(true);
		fuzzyNucleotidePredictionAndEntropy.setSelected(false);
		useSamplingAndAveragingRNAalifoldButton.setSelected(false);
		useSamplingAndAveragingRNAalifoldButton.addChangeListener(this);
		rnaalifoldOptions.setEnabled(false);
		rnaalifoldOptions.setBorder(BorderFactory.createTitledBorder("RNAalifold settings"));
		optionsPanel.add(useSamplingAndAveragingButton, c);
		//optionsPanel.add(fuzzyNucleotidePredictionAndEntropy, c);
		optionsPanel.add(useSamplingAndAveragingRNAalifoldButton, c);
		//optionsPanel.add(consensusEvolutionPrediction, c);
		//cp.setB

		//c.gridx = 0;
		//c.gridy = 0;
		add(optionsPanel, BorderLayout.NORTH);
		
		c.gridx = 0;
		c.gridy = 3;
		c.insets = new Insets(2,2,2,10);
		rnaalifoldOptions.add(new JLabel("RNAalifold path"), c);
		c.insets = new Insets(2,2,2,2);		
		
		c.gridx = 1;
		c.gridy = 3;
		c.gridwidth = 2;
		executableField.setEditable(false);
		rnaalifoldOptions.add(executableField, c);
		//c.gridx = 2;
		//c.gridwidth = 1;
		//pan.add(new JPanel(), c);
		c.gridx = 3;
		c.gridy = 3;
		c.gridwidth = 1;
		c.weightx = 0.5;
		executableButton.addActionListener(this);
		rnaalifoldOptions.add(executableButton, c);
		//pan.add(executablePanel);		
		

		c.gridx = 0;
		c.gridy = 4;
		rnaalifoldOptions.add(new JLabel("Temperature (celsius)"), c);
		temperatureSpinner.addKeyListener(this);
		c.gridx = 1;
		c.gridy = 4;
		rnaalifoldOptions.add(temperatureSpinner, c);
		

		c.gridx = 0;
		c.gridy = 5;
		rnaalifoldOptions.add(new JLabel("Conformation"), c);
		JPanel conformationPanel = new JPanel();
		GridLayout l2 = new GridLayout(1,2);
		l2.setHgap(5);
		l2.setVgap(2);
		conformationPanel.setLayout(l2);
		ButtonGroup conformationGroup = new ButtonGroup();
		conformationGroup.add(linearButton);
		conformationGroup.add(circularButton);
		linearButton.setSelected(true);
		conformationPanel.add(linearButton);
		conformationPanel.add(circularButton);
		c.gridx = 1;
		c.gridy = 5;
		rnaalifoldOptions.add(conformationPanel, c);
		
		c.gridx = 0;
		c.gridy = 6;
		rnaalifoldOptions.add(new JLabel("Covariance term"), c);
		c.gridx = 1;
		c.gridy = 6;
		rnaalifoldOptions.add(covarianceSpinner, c);
		
		c.gridx = 0;
		c.gridy = 7;
		rnaalifoldOptions.add(new JLabel("Non-compatible penalty"), c);
		c.gridx = 1;
		c.gridy = 7;
		rnaalifoldOptions.add(nonCompatibleSpinner, c);
		
		
		c.gridx = 0;
		c.gridy = 8;
		rnaalifoldOptions.add(new JPanel(), c);
		
//		pan.add(new JLabel("Output file:"));
//		pan.add(outFile);
		bigBox.add(rnaalifoldOptions);		
		Box box = Box.createHorizontalBox();
		JButton butt;
		JButton defaultsButton = new JButton("Use defaults");
		defaultsButton.setActionCommand("DEFAULTS");
		defaultsButton.addActionListener(this);
		box.add(defaultsButton);
		box.add(Box.createHorizontalGlue());
		box.add(butt=new JButton("OK"));
		butt.addActionListener(this);
		getRootPane().setDefaultButton(butt);
		box.add(Box.createHorizontalStrut(20));
		box.add(butt=new JButton("Cancel"));
		butt.addActionListener(this);
		bigBox.add(box);
		add(bigBox, SwingConstants.CENTER);
		add(Box.createHorizontalStrut(20), BorderLayout.LINE_START);
		add(Box.createHorizontalStrut(20), BorderLayout.LINE_END);
		//cp.add(Box.createVerticalStrut(15), BorderLayout.PAGE_START);
		add(Box.createVerticalStrut(15), BorderLayout.PAGE_END);
		addKeyListener(this);
		pack();
		loadOptions();
		updateFoldingParametersAndTest();
		setEnabled(rnaalifoldOptions, useSamplingAndAveragingRNAalifoldButton.isSelected());
	}
	
	public void saveOptions()
	{
		if(prefs != null)
		{
			prefs.put("USE_SAMPLING_AND_AVERAGING_PPFOLD", Boolean.toString(useSamplingAndAveragingButton.isSelected()));
			prefs.put("USE_SAMPLING_AND_AVERAGING_RNAALIFOLD", Boolean.toString(useSamplingAndAveragingRNAalifoldButton.isSelected()));
			prefs.put("USE_FUZZY_NUCLEOTIDE", Boolean.toString(fuzzyNucleotidePredictionAndEntropy.isSelected()));
			prefs.put("RNAALIFOLD_EXECUTABLE", executableField.getText());
			prefs.put("RNAALIFOLD_TEMPERATURE", temperatureSpinner.getValue().toString());
			prefs.put("IS_LINEAR", Boolean.toString(linearButton.isSelected()));
			prefs.put("COVARIANCE_VALUE", covarianceSpinner.getValue().toString());
			prefs.put("NON_COMPATABILITY_VALUE", nonCompatibleSpinner.getValue().toString());
		}
		/*try
		{
			BufferedWriter buffer = new BufferedWriter(new FileWriter("rna.options"));
			buffer.write(Boolean.toString(useSamplingAndAveragingButton.isSelected())+"\n");
			buffer.write(Boolean.toString(useSamplingAndAveragingRNAalifoldButton.isSelected())+"\n");
			buffer.write(Boolean.toString(fuzzyNucleotidePredictionAndEntropy.isSelected())+"\n");
			buffer.write(executableField.getText()+"\n");
			buffer.write(((Double)temperatureSpinner.getValue())+"\n");
			buffer.write(linearButton.isSelected()+"\n");
			buffer.write(((Double)covarianceSpinner.getValue())+"\n");
			buffer.write(((Double)nonCompatibleSpinner.getValue())+"\n");
			buffer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}*/
	}
	
	public void loadOptions()
	{
		if(prefs != null)
		{			
			useSamplingAndAveragingButton.setSelected(prefs.getBoolean("USE_SAMPLING_AND_AVERAGING_PPFOLD", true));
			useSamplingAndAveragingRNAalifoldButton.setSelected(prefs.getBoolean("USE_SAMPLING_AND_AVERAGING_RNAALIFOLD", false));
			fuzzyNucleotidePredictionAndEntropy.setSelected(prefs.getBoolean("USE_FUZZY_NUCLEOTIDE", false));
			String defaultExec = "";
			if(System.getProperty("os.name").toLowerCase().contains("windows"))
			{
				File execFile = new File("RNAalifold.exe");
				if(execFile.exists())
				{
					defaultExec = execFile.getAbsolutePath();
				}
			}
			executableField.setText(prefs.get("RNAALIFOLD_EXECUTABLE", defaultExec));
			temperatureSpinner.setValue(prefs.getDouble("RNAALIFOLD_TEMPERATURE", 37.0));
			boolean linear = prefs.getBoolean("IS_LINEAR", true);
			linearButton.setSelected(linear);
			circularButton.setSelected(!linear);
			covarianceSpinner.setValue(prefs.getDouble("COVARIANCE_VALUE", 1.0));
			nonCompatibleSpinner.setValue(prefs.getDouble("NON_COMPATABILITY_VALUE", 1.0));
			//prefs.
		}
		
		/*try
		{
			BufferedReader buffer = new BufferedReader(new FileReader("rna.options"));
			useSamplingAndAveragingButton.setSelected(Boolean.parseBoolean(buffer.readLine()));
			useSamplingAndAveragingRNAalifoldButton.setSelected(Boolean.parseBoolean(buffer.readLine()));
			fuzzyNucleotidePredictionAndEntropy.setSelected(Boolean.parseBoolean(buffer.readLine()));
			executableField.setText(buffer.readLine());
			temperatureSpinner.setValue(Double.parseDouble(buffer.readLine()));
			boolean linear = Boolean.parseBoolean(buffer.readLine());
			linearButton.setSelected(linear);
			circularButton.setSelected(!linear);
			covarianceSpinner.setValue(Double.parseDouble(buffer.readLine()));
			nonCompatibleSpinner.setValue(Double.parseDouble(buffer.readLine()));
			buffer.close();
		}
		catch(IOException ex)
		{

		}*/	
	}
	
	public void useDefaultOptions()
	{
		if(prefs != null)
		{
			try {
				prefs.clear();
			} catch (BackingStoreException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			loadOptions();
		}
	}
	
	public void setEnabled(Container component, boolean enabled)
	{
		Component [] components = ((Container) component).getComponents();
		//component.setEnabled(enabled);
		for(int i = 0 ; i < components.length ;i++)
		{
			components[i].setEnabled(enabled);
			/*if(components[i] instanceof JComponent)
			{				
				((JComponent)components[i]).setOpaque(false);
			}*/
			if(components[i] instanceof Container)
			{				
				setEnabled((Container)components[i], enabled);
			}
		}
	}
	
	
	
	
	void display(Component c) {
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
	@Override
	public void actionPerformed(ActionEvent ev) {
		try{
			if(ev.getSource().equals(executableButton))
			{
				fileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
				int returnVal = fileChooser.showOpenDialog(this);
				if(returnVal == JFileChooser.APPROVE_OPTION)
				{
					File selectedFile = fileChooser.getSelectedFile();
					if(selectedFile.isDirectory())
					{
						File [] files = selectedFile.listFiles();
						for(int i = 0 ; i < files.length ; i++)
						{
							if(files[i].getName().toLowerCase().contains("rnaalifold"))
							{
								executableField.setText(files[i].getAbsolutePath());
								break;
							}
						}
					}
					else
					{
						executableField.setText(selectedFile.getAbsolutePath());
					}
				}				
			}			
			if(ev.getActionCommand().equals("DEFAULTS"))
			{
				useDefaultOptions();
			}
			if(ev.getActionCommand() == "OK") {	

				updateFoldingParametersAndTest();
				this.dispose();
			}
			else
			if(ev.getActionCommand() == "Cancel")
			{
				//setVisible(false);
				this.dispose();
			}
		} catch(NumberFormatException e){
			ErrorMessage.showPane(owner, "Wrong format, "+e.getLocalizedMessage(), false);
		}
	}
	
	public void updateFoldingParametersAndTest()
	{
		boolean sel = useSamplingAndAveragingRNAalifoldButton.isSelected();
		updateFoldingParameters();
		if(sel)
		{
			boolean isWorking = RNAalifold.checkRNAalifold();
			if(sel && !isWorking)
			{
				useSamplingAndAveragingRNAalifoldButton.setSelected(false);
				JOptionPane.showMessageDialog(this,
					    "Disabling RNAalifold folding, the executable does not appear to be working.\nCheck the path or download a newer version of RNAalifold.",
					    "Warning",
					    JOptionPane.WARNING_MESSAGE);					
			}
		}
		updateFoldingParameters();
		saveOptions();
	}
	
	public static void main(String [] args)
	{
		JFrame main = new JFrame();
		main.setSize(new Dimension(900,650));
		main.setLocation(100, 50);
		main.setVisible(true);
		RNASettingsDlg dlg = new RNASettingsDlg(main);
		dlg.display(main);
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

	@Override
	public void stateChanged(ChangeEvent e) {
		if(e.getSource().equals(useSamplingAndAveragingRNAalifoldButton))
		{
			this.setEnabled(rnaalifoldOptions, useSamplingAndAveragingRNAalifoldButton.isSelected());
		}
		
	}
	
	public void updateFoldingParameters()
	{
		
		if(useSamplingAndAveragingButton.isSelected())
		{
			String ppfold = "";
			Postprocess.pluginParameters.setParameter("ppfold", "");
		}
		else
		{
			Postprocess.pluginParameters.removeParameter("ppfold");
		}
		
		if(useSamplingAndAveragingRNAalifoldButton.isSelected())
		{
			RNAalifold.executable = this.executableField.getText();
			String rnaalifold = RNAalifold.executable + " -T " + temperatureSpinner.getValue() +" --cfactor " + covarianceSpinner.getValue() + " --nfactor " + nonCompatibleSpinner.getValue() + " ";
			if(this.circularButton.isSelected())
			{
				rnaalifold += " -circ ";
			}
			Postprocess.pluginParameters.setParameter("rnaalifold", rnaalifold);
		}
		else
		{
			Postprocess.pluginParameters.removeParameter("rnaalifold");
		}
	}
	
}
