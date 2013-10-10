package statalign.ui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Font;
import java.awt.SystemColor;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;

public class RNAPopup extends JDialog implements ActionListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	static String text = "The program has detected that these are RNA sequences. To enable RNA mode, toggle the RNA icon on the toolbar.";
	JButton okButton;
	
	public static void showPane(JFrame owner) {
		JOptionPane.showMessageDialog(owner, printed(text));
	}
	
	public RNAPopup(JFrame owner) {
		super(owner, "RNA Mode", true);
		Container cp = getContentPane();
		cp.setLayout(new BorderLayout(2,1));
		
		
		JTextArea ta = new JTextArea(printed(text));
		ta.setFont(new Font("Monospaced",Font.PLAIN,14));
    	ta.setEnabled(false);
    	ta.setCaretPosition(0);
    	
		cp.add(ta, "North");
		cp.setBackground(SystemColor.controlHighlight);
		
		JPanel panel = new JPanel();
		okButton = new JButton("OK");
		okButton.setActionCommand("OK");
		okButton.addActionListener(this);
		panel.add(okButton);
		
		
		cp.add(panel);
		pack();
		setVisible(true);
		//this.setSize(200, 100);
		
		//setVisible(true);
		
	}
	
	void display(Component c) {
		setLocationRelativeTo(c);
		//System.out.println("Popping up.");
		setSize(600,600);
		this.setResizable(false);
		pack();
		setVisible(true);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		// TODO Auto-generated method stub
		
		if(e.getActionCommand() == "OK") {
			dispose();
		}
		
	}
	
	static String printed(String message){
		String p = "";
		int curpos = 0;
		int nextpos = curpos;
		int prevpos = nextpos;
		while(curpos < message.length()){
		    while((nextpos = nextWhiteSpace(message,prevpos)) < message.length() && nextpos - curpos < 50){
			prevpos = nextpos + 1;
		    }
		    if(nextpos < message.length()){
			if(curpos != prevpos){
			    p += message.substring(curpos,prevpos)+'\n';
			    curpos = prevpos;
			}
			else{
			    p += message.substring(curpos,nextpos)+'\n';
			    curpos = nextpos;
			}
		    }
		    else{
			if(nextpos - curpos < 50){
			    p += message.substring(curpos,nextpos)+'\n';
			    curpos = nextpos+1;		
			}
			else{		    
			    p += message.substring(curpos,prevpos)+'\n'+message.substring(prevpos,nextpos);
			    curpos = nextpos+1;		
			}
		    }
		}
		//	System.out.println(p);
		return p;
	    }

	    static int nextWhiteSpace(String m, int curpos){
		int i = curpos;
		while(i < m.length() && m.charAt(i) != ' ' && m.charAt(i) != '\t' && m.charAt(i) != '\n'){
		    i++;
		}
		return i;
	    }
	
	
}
