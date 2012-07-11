package com.ppfold.main;

import javax.swing.JTextArea;

public class AutoScrollingTextArea extends JTextArea {
	private static final long serialVersionUID = 1L;

	@Override
	public void append(String text) {
		super.append(text);
		this.setCaretPosition(this.getDocument().getLength());  
	}
	
	
}
