package com.ppfold.main;

import java.io.OutputStream;

import javax.swing.JTextArea;
import javax.swing.Scrollable;

public class TextAreaOutputStream extends OutputStream {
	JTextArea textArea;
	
	public TextAreaOutputStream(JTextArea textArea){
		this.textArea = textArea;
	}
	public void flust(){
		textArea.repaint();
	}
	public void write(int b){
		textArea.append(new String(new byte[] {(byte)b}));
	}
	
}
