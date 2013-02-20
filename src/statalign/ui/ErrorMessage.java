package statalign.ui;

import javax.swing.JFrame;
import javax.swing.JOptionPane;


/**
 * A dialog window for showing error messages.
 * All handled error messages that users have to be notified about are shown on the screen by
 * this class.
 */
public class ErrorMessage  {
    
    protected static String except2String(Exception ex) {
			StackTraceElement[] ste = ex.getStackTrace();
			String s = "";
			for(int i = 0; i < ste.length; i++){
				s += ste[i].toString()+"\n";
			}
			return ex.toString()+" in \n"+s;
    }
    
    public static void showPane(JFrame owner, Exception ex, boolean error) {
    	showPane(owner, except2String(ex), error);
    }

    public static void showPane(JFrame owner, String message, String title, boolean error) {
    	JOptionPane.showMessageDialog(owner, message, title == null ? (error ? "Error" : "Message") : title, error ? JOptionPane.ERROR_MESSAGE : JOptionPane.INFORMATION_MESSAGE);
    }
    
    public static void showPane(JFrame owner, String message, boolean error) {
    	showPane(owner, message, null, error);
    }
    
}
