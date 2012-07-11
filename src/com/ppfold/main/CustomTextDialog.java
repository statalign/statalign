package com.ppfold.main;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

public class CustomTextDialog extends JDialog {

	private static final long serialVersionUID = 1L;

	private String text;
	private String title;
	private String message;
	private JTextArea textArea;
	
	
	/**
	 * Generates a dialog window with a custom (not editable) text message, and a text that can be copied to the clipboard.
	 * 
	 * @param title - the title of the window
	 * @param message - the not editable message
	 * @param text - the message that can be copied
	 */
	public CustomTextDialog(String title, String message, String text){
		this.title = title; 
		this.text = text;
		this.message = message+"\n";
		createAndShowDialog();
	}
	
	private void createAndShowDialog(){
		//This builds the dialog window
		setTitle(this.title);
		setLocation(250,40);
		
		Container container = this.getContentPane();
		JPanel contentPane = new JPanel();
		container.add(contentPane);
		contentPane.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		
		contentPane.setLayout(new BoxLayout(contentPane,BoxLayout.PAGE_AXIS));
		
		JTextArea maintext = new JTextArea(message); 
		maintext.setEditable(false);
		maintext.setBackground(contentPane.getBackground());
		contentPane.add(maintext);

		textArea = new JTextArea();
		textArea.setText(this.text);
		textArea.setCaretPosition(0);
		textArea.setEditable(false);
		textArea.setLineWrap(true);
		textArea.setWrapStyleWord(true);

		JScrollPane scrollpane = new JScrollPane(textArea);
		scrollpane.setVerticalScrollBarPolicy(
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		scrollpane.setHorizontalScrollBarPolicy(
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scrollpane.setPreferredSize(new Dimension(100, 100));
		
		contentPane.add(scrollpane);
		
		JButton okButton = new JButton("OK");
		JButton copyButton = new JButton("Copy to clipboard");
		
		JPanel buttonPanel = new JPanel();
		buttonPanel.add(okButton);
		buttonPanel.add(copyButton);
		contentPane.add(buttonPanel);

	    okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	dispose();
            }
        }); 
	    
	    copyButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
        		String selection = textArea.getText();
       		    StringSelection data = new StringSelection(selection);
       		   Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
       		   clipboard.setContents(data, data);
            }
        }); 
	    
		setResizable(false);
		setModal(true);
		pack();
		setVisible(true);
	}

	

}
