package statalign.ui;

import javax.swing.*;
//import javax.swing.filechooser.FileFilter;
//import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
//import java.io.*;


/**
 * A dialog window for showing error messages.
 * All handled error messages that users have to be notified about are shown on the screen by
 * this class.
 */
public class ErrorMessage extends JDialog implements ActionListener {
    
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	Container cp = getContentPane();
    private JButton btOK;

    protected static String except2String(Exception ex) {
			StackTraceElement[] ste = ex.getStackTrace();
			String s = "";
			for(int i = 0; i < ste.length; i++){
				s += ste[i].toString()+"\n";
			}
			return ex.toString()+" in \n"+s;
    }
    
    public ErrorMessage(JFrame owner, Exception ex, boolean error) {
			this(owner, except2String(ex), error);
    }

    public static void showPane(JFrame owner, Exception ex, boolean error) {
    	showPane(owner, except2String(ex), error);
    }

    public static void showPane(JFrame owner, String message, boolean error) {
    	JOptionPane.showMessageDialog(owner, message, error ? "Error" : "Message", error ? JOptionPane.ERROR_MESSAGE : JOptionPane.INFORMATION_MESSAGE);
    }
    
    /**
     * It constructs a dialog window.
     * 
     * @param owner The owner frame. The owner cannot get the focus until the dialog window is closed.
     * @param message The message to be displayed on the dialog window.
     * @param error If true, it is an error message, if false then it is a message.
     */
    public ErrorMessage(JFrame owner, String message, boolean error){
    	super(owner,error ? "Error!!!" : "Message",true);
    	String printedMessage = printed(message);
    	JTextArea ta;
    	setLocation(getParent().getX()+80,getParent().getY()+60);
    	setSize(600,600);
    	setResizable(false);
    	cp.setBackground(SystemColor.controlHighlight);
    	cp.setLayout(new BorderLayout(2,1));
    	cp.add(ta = new JTextArea(printedMessage),"North");
    	ta.setFont(new Font("Monospaced",Font.PLAIN,14));
    	ta.setEnabled(false);
    	ta.setCaretPosition(0);
    	//System.out.println(message);
    	//	cp.add(new JLabel("Can I see this???"));
    	JPanel panel = new JPanel();
    	btOK = new JButton("OK");
    	panel.add(btOK);
    	cp.add(panel);
    	btOK.addActionListener(this);
    	pack();
    	setVisible(true);
    }

    /**
     * This function is for closing the window
     */
    public void actionPerformed(ActionEvent e){
    	dispose();
    }

    String printed(String message){
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

    int nextWhiteSpace(String m, int curpos){
	int i = curpos;
	while(i < m.length() && m.charAt(i) != ' ' && m.charAt(i) != '\t' && m.charAt(i) != '\n'){
	    i++;
	}
	return i;
    }

}
