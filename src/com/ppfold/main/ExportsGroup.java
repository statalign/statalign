package com.ppfold.main;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

public class ExportsGroup extends JPanel implements ActionListener {

	private static final long serialVersionUID = -6196506765119020438L;
	private JRadioButton onlyCTButton;
	private JRadioButton defaultButton;
	private JRadioButton exportAllButton;
	private JTextField exportNameField;
	
	public ExportsGroup(){
		this.setLayout(new BoxLayout(this,BoxLayout.LINE_AXIS));
		
		onlyCTButton = new JRadioButton("Only .ct");
		defaultButton = new JRadioButton("5 files");
		exportAllButton = new JRadioButton("Export all");
		
		onlyCTButton.setToolTipText("Export only the .ct file");
		defaultButton.setToolTipText("Export .ct, .st, .lseq and .seq files and .newick if tree is not given");
		exportAllButton.setToolTipText("As the default, plus basepairing matrices (might generate very large files!)");
	
		onlyCTButton.addActionListener(this);
		defaultButton.addActionListener(this);
		exportAllButton.addActionListener(this);
		
		ButtonGroup group = new ButtonGroup();
		group.add(onlyCTButton);
		group.add(defaultButton);
		group.add(exportAllButton);
			
		JLabel label = new JLabel("Exports:");
		add(label);
		add(onlyCTButton);
		add(defaultButton);
		add(exportAllButton);
		onlyCTButton.setSelected(true);
		PPfoldMain.onlyCT = true;
		
		JLabel nextlabel = new JLabel("     File prefix: ");
		exportNameField = new JTextField("");
		exportNameField.setToolTipText("Leave empty to match alignment name.");
		exportNameField.setPreferredSize(new Dimension(100,20));
		add(nextlabel);
		add(exportNameField);
		
		exportNameField.getDocument().addDocumentListener(new DocumentListener(){
			public void changedUpdate(DocumentEvent arg0) {updateExportFileHandle();}
			public void insertUpdate(DocumentEvent arg0) {updateExportFileHandle();}
			public void removeUpdate(DocumentEvent arg0) {updateExportFileHandle();}
		});
		
	}

	public void updateExportNameField(String name){
		exportNameField.setText(name);
	}
	
	public void updateExportFileHandle(){
		if(!exportNameField.getText().isEmpty()){
			PPfoldMain.exportfilehandle = exportNameField.getText();
			PPfoldMain.specialname = true;
		}
		else{		
			PPfoldMain.specialname = false;
			PPfoldMain.exportfilehandle = null;
		}
	}
	
	public void actionPerformed(ActionEvent arg0) {
		if(onlyCTButton.isSelected()){
			PPfoldMain.onlyCT = true;
			PPfoldMain.exportson = false;
		}
		else if(defaultButton.isSelected()){
			PPfoldMain.onlyCT = false;
			PPfoldMain.exportson = false;
		}
		else if(exportAllButton.isSelected()){
			PPfoldMain.onlyCT = false;
			PPfoldMain.exportson = true;
		}
		else{
			defaultButton.setSelected(true);
			PPfoldMain.onlyCT = false;
			PPfoldMain.exportson = false;
		}
	}
	@Override 
	public void setEnabled(boolean value){
		onlyCTButton.setEnabled(value);
		defaultButton.setEnabled(value);
		exportAllButton.setEnabled(value);
		exportNameField.setEnabled(value);
	}
	
	public String getExportPrefix(){
		return exportNameField.getText();
	}

	public void setExportPrefix(String exportfilehandle) {
		exportNameField.setText(exportfilehandle);
		updateExportFileHandle();
	}
}
