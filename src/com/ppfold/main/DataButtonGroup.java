package com.ppfold.main;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

public class DataButtonGroup extends JPanel implements Updatable {
	  /**
	 * 
	 */
	private static final long serialVersionUID = 8735013149334446736L;
	// members:
	  private JLabel title; 
	  private JButton addButton;
	  private JButton editButton;
	  private JButton removeButton;
	  private JComboBox chooser;
	  JTextArea dataInfoTextArea ;
	  
	  public static int FILECHOOSER = 0;
	  public static int FOLDERCHOOSER = 1;
	  
	  // constructors:
	  public DataButtonGroup() {

		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		//this.setSize(100,250);
		//this.setPreferredSize(new Dimension(100,250));
		
	    // create components
		JPanel comboPanel = new JPanel();
		comboPanel.setLayout(new BoxLayout(comboPanel, BoxLayout.LINE_AXIS));
		//comboPanel.setSize(new Dimension(100, 20));
		//comboPanel.setPreferredSize(new Dimension(100, 20));
		comboPanel.setAlignmentX(Component.CENTER_ALIGNMENT);
		title = new JLabel("Data and constraints: ");
		//title.setPreferredSize(new Dimension(130, 10));
		title.setAlignmentX(Component.CENTER_ALIGNMENT);
		title.setHorizontalAlignment(JLabel.CENTER);
		comboPanel.add(title);
		chooser = new JComboBox();
		//chooser.setPreferredSize(new Dimension(200, 15));
		chooser.setEnabled(false);
		chooser.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	int itemnr = chooser.getItemCount();
            	if(itemnr>0){
            		int nr = chooser.getSelectedIndex();
            		DataInfo datainfo = PPfoldMain.datainfo.get(nr);
            		String text = "";
            		if(datainfo.getType()==0){
            			text = "File: " + datainfo.getFileName() + "\n";
            			text += "Distribution: " + ((datainfo.getDistFileName())==
            					PPfoldMain.defaultDataDistfile?"Default\n":(datainfo.getDistFileName() + "\n"));
            			text += "Sequence name: " + datainfo.getSequenceName();
            		}
            		else if(datainfo.getType()==2){
            			//Constraint data
            			if(datainfo.getFileName()!=null){
            				text = "File: " + datainfo.getFileName() + "\n";
            			}
            			if(datainfo.getSequenceName()!=null){
                			text += "Sequence name: " + datainfo.getSequenceName() + "\n";            			
                			}
            			if(datainfo.getContactDistance()>-1){
            				text += "Max. contact distance: " + datainfo.getContactDistance() + "\n";
            			}

            		}
            		else{
            			//Constraint data
            			text = "File: " + datainfo.getFileName() + "\n";
            			text += "Sequence name: " + datainfo.getSequenceName() + "\n";
            		}
            		dataInfoTextArea.setText(text);
            		
            	}
            	else{
            		dataInfoTextArea.setText("<data not selected>");
            	}
            }
        });
	    
		chooser.setAlignmentX(Component.CENTER_ALIGNMENT);
		comboPanel.add(chooser);
		chooser.setPreferredSize(chooser.getSize()); //To avoid resizing when items are added. 
		
		addButton = new JButton("Add");
		addButton.setToolTipText("Add data");
	    //addButton.setPreferredSize(new Dimension(90,30));

	    final DataButtonGroup bg = this;
	    addButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	DataChooserDialog dia = new DataChooserDialog(bg);
            	if(chooser.getItemCount()>0){
            		removeButton.setEnabled(true);
            		chooser.setEnabled(true);
            		editButton.setEnabled(true);
            	}
            }
        });
	    
	    editButton = new JButton("Edit");
		editButton.setToolTipText("Edit data");
	    //addButton.setPreferredSize(new Dimension(90,30));
		editButton.setEnabled(false);
	    editButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	DataChooserDialog dia = new DataChooserDialog(PPfoldMain.datainfo.get(chooser.getSelectedIndex()), bg);
            	if(chooser.getItemCount()>0){
            		removeButton.setEnabled(true);
            		chooser.setEnabled(true);
            	}
            }
        });
	    
	    
	    
		removeButton = new JButton("Remove");
		removeButton.setToolTipText("Remove data");
		removeButton.setEnabled(false);
	   // removeButton.setPreferredSize(new Dimension(80,30));
	    removeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	int selected = chooser.getSelectedIndex();
            	PPfoldMain.datainfo.remove(selected);
            	updateGUI();
            	if(chooser.getItemCount()==0){
            		removeButton.setEnabled(false);
            		chooser.setEnabled(false);
            		editButton.setEnabled(false);
            		PPfoldMain.auxdata=false;
            	}
            }
        });
	    

		comboPanel.add(Box.createRigidArea(new Dimension(5, 0)));
		comboPanel.add(addButton);
		comboPanel.add(Box.createRigidArea(new Dimension(5, 0)));
		comboPanel.add(editButton);
		comboPanel.add(Box.createRigidArea(new Dimension(5, 0)));
		comboPanel.add(removeButton);
		
		JPanel infoPanel = new JPanel();
		infoPanel.setLayout(new BoxLayout(infoPanel, BoxLayout.PAGE_AXIS));
		infoPanel.setAlignmentX(Component.CENTER_ALIGNMENT);
		
		dataInfoTextArea = new JTextArea();
		infoPanel.add(dataInfoTextArea);
		dataInfoTextArea.setText("<data not selected>");
		dataInfoTextArea.setPreferredSize(new Dimension(400, 80));
		dataInfoTextArea.setLineWrap(true);
		float[] somecolor = new float[3];
		somecolor = Color.RGBtoHSB(152,251,152	, null);
		//somecolor = Color.RGBtoHSB(255,255,240, null);
		dataInfoTextArea.setBackground(Color.getHSBColor(somecolor[0], somecolor[1], somecolor[2]));
		//infoPanel.setBackground(Color.BLUE);
		this.add(comboPanel);
		this.add(Box.createRigidArea(new Dimension(0, 5)));
		this.add(infoPanel);
		this.add(Box.createRigidArea(new Dimension(0, 5)));
		this.setVisible(true);
	  }
	  
	  public boolean checkData(){
		  System.out.println("Checking data...");
		  return true;
	  }

	public void updateGUI() {
		int n = PPfoldMain.datainfo.size();
		String[] datalabels = new String[n];
		for(int i = 0; i<n; i++){
			String datatype = "data";
			DataInfo datainfo = PPfoldMain.datainfo.get(i);
			if(datainfo.getType()==0){
				datatype = "(Probing data)";
			}
			else if(datainfo.getType()==1){
				datatype = "(Probability data)";
			}
			else{
				datatype = "(Constraint data)";
			}
			datalabels[i] = new String(PPfoldMain.datainfo.get(i).getiD() + " " + datatype);
		}
	    final DefaultComboBoxModel model = new DefaultComboBoxModel(datalabels);
		chooser.setModel(model);
		chooser.setSelectedIndex(n-1);
		chooser.updateUI();
	}
	
	  @Override 
	  public void setEnabled(boolean value){
		  chooser.setEnabled(value);
		  addButton.setEnabled(value);
		  removeButton.setEnabled(value);
		  //Just to make sure the remove button doesn't get re-enabled
		  if(chooser.getItemCount()==0){
			  removeButton.setEnabled(false);
			  chooser.setEnabled(false);
		  }
		  dataInfoTextArea.setEnabled(value);
	  }

	public void updateModel() {}

	public int getDataNumber() {
		return chooser.getItemCount();
	}

	} 