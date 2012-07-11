package com.ppfold.main;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.File;
import java.net.URL;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.ParseException;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
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
import javax.swing.JFormattedTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.JFormattedTextField.*;
import javax.swing.text.NumberFormatter;

public class DataChooserDialog extends JDialog {

	private static final long serialVersionUID = 8528742785435875483L;

	private JButton okButton;
	private JButton cancelButton;
	private JTextField fileName;
	private JTextField fileName2;
	private JTextField fileName3;
	private JTextField identifier;
	private JTextField seqName;
	private JTextField distName;
	private JTextField contactDistance;
	private JRadioButton barDataButton;
	private JRadioButton constraintDataButton;
	private JRadioButton exactMappingDataButton; 
	private JPanel barDataPanel;
	private JPanel constraintPanel;
	private JPanel advancedPanel;
	private JCheckBox checkBox;
	private JCheckBox checkBox2; 
	private JButton addDistButton; 
	
	private static boolean editing = false;
	private DataInfo privateData  = null; 
	
	public DataChooserDialog(DataButtonGroup parent){
		editing = false; 
		createDialog(parent);
		showDialog(parent);
	}
	
	public DataChooserDialog(DataInfo datainfo, DataButtonGroup parent)
	{
		editing = true; 
		privateData = datainfo; 
		createDialog(parent);
		//Populate with datainfo details
		switch(datainfo.getType()){
		case 0:
			barDataButton.setSelected(true);
			updateGUI();
			fileName.setText(datainfo.getFileName());
			if(datainfo.getDistFileName().equals(PPfoldMain.defaultDataDistfile)){
				checkBox.setSelected(true);
				distName.setEnabled(false);
				addDistButton.setEnabled(false);
			}
			else{
				checkBox.setSelected(false);
				distName.setText(datainfo.getDistFileName());
				distName.setEnabled(true);
			}
			break;
		case 1:
			exactMappingDataButton.setSelected(true);
			fileName3.setText(datainfo.getFileName());
			updateGUI();
			break;
		case 2:
			constraintDataButton.setSelected(true);
			updateGUI();
			if(datainfo.getFileName()!=null){
				fileName2.setText(datainfo.getFileName());
			}
			if(datainfo.getContactDistance()>0){
				String cdString = "" + datainfo.getContactDistance();
				contactDistance.setText(cdString);
				checkBox2.setSelected(true);
			}
			break;
		default:
			System.err.println("Invalid datatype, ignoring...");
			break;
		}
		if(datainfo.getSequenceName()!=null){
			seqName.setText(datainfo.getSequenceName());
		}
		if(datainfo.getiD()!=null){
			identifier.setText(datainfo.getiD());
		}

		showDialog(parent);
		
	}
	
	private void createDialog(final DataButtonGroup parent){
		//This builds the dialog window
		
		setTitle("Data setup");
		setLocation(250,40);
		setResizable(false);
		
		Container container = this.getContentPane();
		JPanel contentPane = new JPanel();
		container.add(contentPane);
		contentPane.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		contentPane.setLayout(new BoxLayout(contentPane,BoxLayout.PAGE_AXIS));
		
		//This panel contains components to allow choosing between different types of data. 
		JPanel radioButtonPanel = new JPanel(); 
		radioButtonPanel.setLayout(new BoxLayout(radioButtonPanel, BoxLayout.LINE_AXIS));
		ButtonGroup dataButtonGroup = new ButtonGroup();
		barDataButton = new JRadioButton("Probing data");
		constraintDataButton = new JRadioButton("Hard constraints");
		exactMappingDataButton = new JRadioButton("Exact probabilities (advanced)");
		dataButtonGroup.add(barDataButton);
		dataButtonGroup.add(constraintDataButton);
		dataButtonGroup.add(exactMappingDataButton);
		radioButtonPanel.add(barDataButton);
		radioButtonPanel.add(constraintDataButton);
		radioButtonPanel.add(exactMappingDataButton);	
		barDataButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				updateGUI();
			}
		});
		constraintDataButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				updateGUI();
			}
		});
		exactMappingDataButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				updateGUI();
			}
		});
		
		
		
		TitledBorder border1 = BorderFactory.createTitledBorder(
				BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Type of data");
		border1.setTitleJustification(TitledBorder.LEFT);
		radioButtonPanel.setBorder(border1);		

		//This panel contains components to select bar data
		barDataPanel = new JPanel();
		barDataPanel.setLayout(new BoxLayout(barDataPanel,BoxLayout.PAGE_AXIS));
		
		JPanel filePanel = new JPanel();
		filePanel.setLayout(new BoxLayout(filePanel,BoxLayout.LINE_AXIS));
		JLabel filenameLabel = new JLabel("Data file:");
		filenameLabel.setPreferredSize(new Dimension(100,10));
		filePanel.add(filenameLabel);
		fileName = new JTextField("<No file selected>");
		fileName.setPreferredSize(new Dimension(121,19));
		fileName.setToolTipText("This file contains the measured experimental data for one sequence of the alignment.");
		fileName.setEditable(false);
		filePanel.add(fileName);
		JButton addFileButton = new JButton("Browse"); 
		addFileButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	final JFileChooser fc = new JFileChooser();
    			if(PPfoldGUIMainWindow.directory!=null){
    				File theDirectory = new File(PPfoldGUIMainWindow.directory);
    			fc.setCurrentDirectory(theDirectory);
    			}
    			int returnVal = fc.showOpenDialog(null);
    	        if (returnVal == JFileChooser.APPROVE_OPTION) {
    	            File file = fc.getSelectedFile();
    	            fileName.setText(file.getAbsolutePath());
            		//PPfoldGUIMainWindow.directory = file.getParent(); ;
            		//System.out.println(PPfoldGUIMainWindow.directory);
    	        }  //otherwise user cancelled, do nothing }     			
    		}            
        }); 
		filePanel.add(addFileButton);
		
		JPanel distPanel = new JPanel();
		distPanel.setLayout(new BoxLayout(distPanel,BoxLayout.LINE_AXIS));
		JLabel distLabel = new JLabel("Distribution file:");
		distLabel.setPreferredSize(new Dimension(100,10));
		distPanel.add(distLabel);
		distName = new JTextField("<No file selected>");
		distName.setEditable(false);
		distName.setPreferredSize(new Dimension(121,19));
		distName.setToolTipText("This file contains the distribution data. See PPfold website for more information.");
		distPanel.add(distName);
		addDistButton = new JButton("Browse"); 
		addDistButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	final JFileChooser fc = new JFileChooser();
    			if(PPfoldGUIMainWindow.directory!=null){
    				File theDirectory = new File(PPfoldGUIMainWindow.directory);
    			fc.setCurrentDirectory(theDirectory);
    			}
    			int returnVal = fc.showOpenDialog(null);
    	        if (returnVal == JFileChooser.APPROVE_OPTION) {
    	            File file = fc.getSelectedFile();
    	            distName.setText(file.getAbsolutePath());
            		//PPfoldGUIMainWindow.directory = file.getParent(); ;
            		//System.out.println(PPfoldGUIMainWindow.directory);
    	        }  //otherwise user cancelled, do nothing }     			
    		}            
        }); 
		checkBox = new JCheckBox("Use default");
		checkBox.setSelected(true);
		checkBox.setToolTipText("The default is the empirical distribution of SHAPE values from E.coli 16S and 23S rRNA.");
		checkBox.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent arg0) {
				if(checkBox.isSelected()){
					addDistButton.setEnabled(false);
					distName.setEnabled(false);
				}
				else{
					addDistButton.setEnabled(true);
					distName.setEnabled(true);
				}
			}
		});
		
		distPanel.add(checkBox);
		distPanel.add(addDistButton);

		barDataPanel.add(filePanel);
		barDataPanel.add(distPanel);
		TitledBorder border2 = BorderFactory.createTitledBorder(
				BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Probing data");
		border2.setTitleJustification(TitledBorder.LEFT);
		//border2.setTitlePosition(TitledBorder.BELOW_TOP);
		barDataPanel.setBorder(border2);	
		
		//This panel contains components to select hard constraint data
		constraintPanel = new JPanel();
		constraintPanel.setLayout(new BoxLayout(constraintPanel,BoxLayout.LINE_AXIS));
		JPanel filePanel2 = new JPanel();
		filePanel2.setLayout(new BoxLayout(filePanel2,BoxLayout.LINE_AXIS));
		JLabel filenameLabel2 = new JLabel("Data file:");
		filenameLabel2.setPreferredSize(new Dimension(100,10));
		filePanel2.add(filenameLabel2);
		fileName2 = new JTextField("<No file selected>");
		fileName2.setPreferredSize(new Dimension(121,19));
		fileName2.setToolTipText("This file contains the constraints in mfold format. See documentation for details");
		fileName2.setEditable(false);
		filePanel2.add(fileName2);
		JButton addFileButton2 = new JButton("Browse"); 
		addFileButton2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	final JFileChooser fc = new JFileChooser();
    			if(PPfoldGUIMainWindow.directory!=null){
    				File theDirectory = new File(PPfoldGUIMainWindow.directory);
    			fc.setCurrentDirectory(theDirectory);
    			}
    			int returnVal = fc.showOpenDialog(null);
    	        if (returnVal == JFileChooser.APPROVE_OPTION) {
    	            File file = fc.getSelectedFile();
    	            fileName2.setText(file.getAbsolutePath());
            		//PPfoldGUIMainWindow.directory = file.getParent(); ;
            		//System.out.println(PPfoldGUIMainWindow.directory);
    	        }  //otherwise user cancelled, do nothing }     			
    		}            
        }); 
		filePanel2.add(addFileButton2);
		constraintPanel.add(filePanel2);
		
		JPanel cdPanel = new JPanel();
		cdPanel.setLayout(new BoxLayout(cdPanel,BoxLayout.LINE_AXIS));
		checkBox2 = new JCheckBox("Max. contact distance: "); 
		checkBox2.setSelected(false);
		cdPanel.add(checkBox2);
		contactDistance = new JTextField();
		contactDistance.addKeyListener(new KeyListener(){
			public void keyPressed(KeyEvent arg0) {}
			public void keyReleased(KeyEvent arg0) {}
			public void keyTyped(KeyEvent arg0) {
				if(!Character.isDigit(arg0.getKeyChar()) && arg0.getKeyChar() != 8 && arg0.getKeyChar() != 127) {
					arg0.consume();
				}
			}
			
		});
		checkBox2.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent arg0) {
				if(checkBox2.isSelected()){
					contactDistance.setEnabled(true);
				}
				else{
					contactDistance.setEnabled(false);
				}
			}
		});
		checkBox2.setToolTipText("Specify the maximum distance between basepairs. (Might not be exactly fulfilled in highly gapped regions.)");
		contactDistance.setToolTipText("Specify the maximum distance between basepairs. (Might not be exactly fulfilled in highly gapped regions.)");

		contactDistance.setEnabled(false);
		cdPanel.add(contactDistance);
		constraintPanel.add(cdPanel);
		
		
		TitledBorder border3 = BorderFactory.createTitledBorder(
				BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Hard constraints");
		border3.setTitleJustification(TitledBorder.LEFT);
		//border2.setTitlePosition(TitledBorder.BELOW_TOP);
		constraintPanel.setBorder(border3);	

		//This panel contains components to select probability data
		advancedPanel = new JPanel();
		advancedPanel.setLayout(new BoxLayout(advancedPanel,BoxLayout.LINE_AXIS));
		JPanel filePanel3 = new JPanel();
		filePanel3.setLayout(new BoxLayout(filePanel3,BoxLayout.LINE_AXIS));
		JLabel filenameLabel3 = new JLabel("Data file:");
		filenameLabel3.setPreferredSize(new Dimension(100,10));
		filePanel3.add(filenameLabel3);
		fileName3 = new JTextField("<No file selected>");
		fileName3.setToolTipText("This file contains the probability constraints.");
		fileName3.setEditable(false);
		fileName3.setPreferredSize(new Dimension(121,19));
		filePanel3.add(fileName3);
		JButton addFileButton3 = new JButton("Browse"); 
		addFileButton3.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				final JFileChooser fc = new JFileChooser();
				if(PPfoldGUIMainWindow.directory!=null){
					File theDirectory = new File(PPfoldGUIMainWindow.directory);
					fc.setCurrentDirectory(theDirectory);
				}
				int returnVal = fc.showOpenDialog(null);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					fileName3.setText(file.getAbsolutePath());
				}  //otherwise user cancelled, do nothing }     			
			}            
		}); 
		filePanel3.add(addFileButton3);
		advancedPanel.add(filePanel3);
		TitledBorder border4 = BorderFactory.createTitledBorder(
				BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Exact probabilities");
		border4.setTitleJustification(TitledBorder.LEFT);
		//border2.setTitlePosition(TitledBorder.BELOW_TOP);
		advancedPanel.setBorder(border4);	
		
		JPanel seqnamePanel = new JPanel();
		JLabel nameLabel = new JLabel("Enter sequence name:");
		//nameLabel.setPreferredSize(new Dimension(120,20));
		seqName = new JTextField(""); 
		seqName.setToolTipText("Name of the sequence for which the data were obtained. Must match the name of a sequence in the alignment.");
		seqName.setPreferredSize(new Dimension(100,20));
		//nameLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
		//seqName.setAlignmentX(Component.LEFT_ALIGNMENT);
		seqnamePanel.add(nameLabel);
		seqnamePanel.add(seqName);

		JLabel idLabel = new JLabel("Enter identifier name (optional):");
		//nameLabel.setPreferredSize(new Dimension(120,20));
		identifier = new JTextField("Data " + (parent.getDataNumber()+1)); 
		identifier.setPreferredSize(new Dimension(100,20));
		identifier.setToolTipText("The name you want to call this dataset. (For your use only, not used by the algorithm.)");
		//nameLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
		//seqName.setAlignmentX(Component.LEFT_ALIGNMENT);
		seqnamePanel.add(idLabel);
		seqnamePanel.add(identifier);		

		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new BoxLayout(buttonPanel,BoxLayout.LINE_AXIS));
		okButton = new JButton("OK");
		buttonPanel.add(okButton);
		okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	DataInfo datainfo = new DataInfo();
            	if(barDataButton.isSelected()){
            		if(!fileName.getText().startsWith("<No file selected>")){
            			datainfo.setFileName(fileName.getText());
            		}
            		else{
            			JOptionPane.showMessageDialog(null,
            					"File name must be given!",
            					"Missing input",
            					JOptionPane.ERROR_MESSAGE);
            			return;
            		}

            		if(checkBox.isSelected()){
            			datainfo.setDistFileName(PPfoldMain.defaultDataDistfile);
            		}
            		else{
            			if(!distName.getText().startsWith("<No file selected>")){
            				datainfo.setDistFileName(distName.getText());
            			}
            			else{
            				JOptionPane.showMessageDialog(null,
            						"Distribution must be given!",
            						"Missing input",
            						JOptionPane.ERROR_MESSAGE);
            				return;
            			}
            		}

            		if(!seqName.getText().trim().isEmpty()){
            			datainfo.setSequenceName(seqName.getText());
            		}
            		else{
            			JOptionPane.showMessageDialog(null,
            					"Sequence name must be given!\n(It must match the name in the alignment)",
            					"Missing input",
            					JOptionPane.ERROR_MESSAGE);
            			return;
            		}
            		datainfo.setiD(identifier.getText());
            		datainfo.setType(0);
    				datainfo.setContactDistance(-1);
            		PPfoldMain.auxdata=true;
            	}
            	else if(constraintDataButton.isSelected()){
            		if(fileName2.getText().startsWith("<No file selected>") && 
            				(!checkBox2.isSelected() || (checkBox2.isSelected() && 
            						contactDistance.getText() != null && Integer.valueOf(contactDistance.getText())<=0))
            						){
            			JOptionPane.showMessageDialog(null,
            					"File name or positive contact distance must be given!",
            					"Missing input",
            					JOptionPane.ERROR_MESSAGE);
            			return;
            		}
            		if(!fileName2.getText().startsWith("<No file selected>") && seqName.getText().trim().isEmpty()){
            			JOptionPane.showMessageDialog(null,
            					"Sequence name must be given if constraint file is selected! \n(It must match the name in the alignment)",
            					"Missing input",
            					JOptionPane.ERROR_MESSAGE);
            			return;
            		}
            		if(checkBox2.isSelected() && (contactDistance.getText()==null || contactDistance.getText().trim().equals(""))){
            			JOptionPane.showMessageDialog(null,
            					"Contact distance is missing!",
            					"Missing input",
            					JOptionPane.ERROR_MESSAGE);
            			return;
            		}
            		
            		if(!fileName2.getText().startsWith("<No file selected>")){
            			datainfo.setFileName(fileName2.getText());
            		}

            		
            		if(checkBox2.isSelected()&&contactDistance.getText()!=null){
            			//Contact distance specified
        				datainfo.setContactDistance(Integer.valueOf(contactDistance.getText()));
            		}
            		else{
            			//No contact distance specified
        				datainfo.setContactDistance(-1);
            		}
            		if(seqName.getText()!=null && !seqName.getText().trim().isEmpty()){
            			datainfo.setSequenceName(seqName.getText());
            		}

            		datainfo.setiD(identifier.getText());
            		datainfo.setType(2);
            		PPfoldMain.auxdata=true;
            	}
            	else{
            		if(!fileName3.getText().startsWith("<No file selected>")){
            			datainfo.setFileName(fileName3.getText());
            		}
            		else{
            			JOptionPane.showMessageDialog(null,
            					"File name must be given!",
            					"Missing input",
            					JOptionPane.ERROR_MESSAGE);
            			return;
            		}
            		if(seqName.getText()!=null && !seqName.getText().trim().isEmpty()){
            			datainfo.setSequenceName(seqName.getText());
            		}
            		else{
            			JOptionPane.showMessageDialog(null,
            					"Sequence name must be given!\n(It must match the name in the alignment)",
            					"Missing input",
            					JOptionPane.ERROR_MESSAGE);
            			return;
            		}
            		datainfo.setiD(identifier.getText());
            		datainfo.setType(1);
            		PPfoldMain.auxdata=true;
    				datainfo.setContactDistance(-1);

            	}
            	PPfoldMain.datainfo.add(datainfo);
        		if(editing){
        			PPfoldMain.datainfo.remove(privateData);
        		}
            	parent.updateGUI();
        		dispose();
            	
            }

        }); 
		
		
	
		buttonPanel.add(new JLabel(" "));
		cancelButton = new JButton("Cancel");
		buttonPanel.add(cancelButton);
		
		cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	dispose();
            }
        }); 
		
		radioButtonPanel.setAlignmentX(Component.CENTER_ALIGNMENT);
		barDataPanel.setAlignmentX(Component.CENTER_ALIGNMENT);
		constraintPanel.setAlignmentX(Component.CENTER_ALIGNMENT);
		advancedPanel.setAlignmentX(Component.CENTER_ALIGNMENT);
		seqnamePanel.setAlignmentX(Component.CENTER_ALIGNMENT);
		buttonPanel.setAlignmentX(Component.CENTER_ALIGNMENT);
		
		contentPane.add(radioButtonPanel);
		contentPane.add(barDataPanel);
		contentPane.add(constraintPanel);
		contentPane.add(advancedPanel);
		contentPane.add(seqnamePanel);
		contentPane.add(buttonPanel);
		
		barDataButton.setSelected(true);
		enablePanel(barDataPanel);
		disablePanel(constraintPanel);
		disablePanel(advancedPanel);
		addDistButton.setEnabled(false);
		distName.setEnabled(false);

	}
	
	private void showDialog(JPanel panel){
		setModal(true);
		pack();
		if(editing){
			switch(privateData.getType()){
			case 0: 
				barDataButton.requestFocusInWindow(); 
				break; 
			case 1:  
				exactMappingDataButton.requestFocusInWindow(); 
				break;
			case 2: 
				constraintDataButton.requestFocusInWindow(); 
				break; 
			default: break; 
			}
		}
		setVisible(true);	
	}
	
	private void disablePanel(JPanel panel){
		if(panel!=null){
			panel.setEnabled(false);
			for(Component comp:panel.getComponents()){
				if(comp instanceof JPanel){
					disablePanel((JPanel) comp);
				}
				comp.setEnabled(false);
			}}
		else{
			System.out.println("Panel was null!");
		}
	}

	private void enablePanel(JPanel panel){
		if(panel!=null){
			panel.setEnabled(true);
			for(Component comp:panel.getComponents()){
				if(comp instanceof JPanel){
					enablePanel((JPanel) comp);
				}
				comp.setEnabled(true);
			}
		}
		else{
			System.err.println("Panel was null!");
		}
	}
	
	private void updateGUI(){
		if(barDataButton.isSelected()){
			enablePanel(barDataPanel);
			//if(checkBox.isSelected()){
			//	addDistButton.setEnabled(false);
			//	distName.setEnabled(false);
			//}
			disablePanel(constraintPanel);
			disablePanel(advancedPanel);
		}
		else if(constraintDataButton.isSelected()){
			disablePanel(barDataPanel);
			enablePanel(constraintPanel);
			disablePanel(advancedPanel);			
		}
		else{
			disablePanel(barDataPanel);
			disablePanel(constraintPanel);
			enablePanel(advancedPanel);
		}
	}
	
}
