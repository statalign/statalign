/**
 * 
 */
package statalign.ui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

/**
 * This is a help window for displaying html pages.
 * 
 * @author miklos
 *
 */
public class HelpWindow extends JDialog implements ActionListener{
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private final int WIDTH = 600;
    private final int HEIGHT = 400;
    private JEditorPane editorpane;
    private URL homeURL;
    private ArrayList<URL> visitedPages;
    private int visitedIndex;
    private JButton btnback;
    private JButton btnforward;
    
/**
 * This constructs a dialog window displaying html messages.
 * @param owner The owner frame of the dialog window. The owner cannot get back the focus until the
 * dialog window is closed.
 * @param title The title of the window.
 * @param hlpURL The home url that is displayed when the window is displayed.
 * @param buttons It true, navigation buttons are shown.
 */
public HelpWindow(JFrame owner, String title, URL hlpURL, boolean buttons) {
    super(owner,title, true);
    homeURL = hlpURL;  
    editorpane = new JEditorPane();
    editorpane.setEditable(false);
    try {
        editorpane.setPage(homeURL);
    } catch (Exception ex) {
        ErrorMessage.showPane(null,"name: "+homeURL+"\n"+ex.getMessage(),true);
    }
    visitedPages = new ArrayList<URL>();
    visitedIndex = 0;
    visitedPages.add(homeURL);
    //anonymous inner listener
    editorpane.addHyperlinkListener(new HyperlinkListener() {
        @Override
		public void hyperlinkUpdate(HyperlinkEvent ev) {
            try {
                if (ev.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                	for(int i = visitedIndex + 1; i < visitedPages.size(); ){
                		visitedPages.remove(i);
                		//System.out.println("Removing index "+i+" size: "+visitedPages.size());
                	}
                	visitedPages.add(ev.getURL());
                	btnback.setEnabled(true);
                	btnforward.setEnabled(false);
                	visitedIndex++;
                  	//System.out.println("visited index: "+visitedIndex+" size of the array: "+visitedPages.size());
                    editorpane.setPage(ev.getURL());
                   // new ErrorMessage(null,"the name of the URL is: "+ev.getURL().getFile(),false);
                }
            } catch (IOException ex) {
                //put message in window
            	ErrorMessage.showPane(null,"name: "+homeURL+"\n"+ex.getMessage(),true);
            }
        }
    });
    getContentPane().add(new JScrollPane(editorpane));
    if(buttons){
    	addButtons();
    }
    // no need for listener just dispose
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    // dynamically set location
    calculateLocation();
    setVisible(true);
    // end constructor
}
/**
 * An Actionlistener is implemented, so must implement this method
 *
 */
@Override
public void actionPerformed(ActionEvent e) {
    String strAction = e.getActionCommand();
   // URL tempURL;
    try {
        if (strAction == "< Back") {
            //       tempURL = editorpane.getPage();
        	btnforward.setEnabled(true);
        	visitedIndex--;
          	//System.out.println("visited index: "+visitedIndex+" size of the array: "+visitedPages.size());
          	editorpane.setPage(visitedPages.get(visitedIndex));
        	if(visitedIndex == 0){
        		btnback.setEnabled(false);
        	}
        }
        if (strAction == "Forward >") {
            //       tempURL = editorpane.getPage();
        	btnback.setEnabled(true);
        	visitedIndex++;
          	//System.out.println("visited index: "+visitedIndex+" size of the array: "+visitedPages.size());
          	editorpane.setPage(visitedPages.get(visitedIndex));
        	if(visitedIndex == visitedPages.size() - 1){
        		btnforward.setEnabled(false);
        	}
        }
       if (strAction == "Home") {
     //       tempURL = editorpane.getPage();
         	for(int i = visitedIndex + 1; i < visitedPages.size(); i++){
        		visitedPages.remove(i);
        	}
        	visitedPages.add(homeURL);
        	btnback.setEnabled(true);
        	btnforward.setEnabled(false);
        	visitedIndex++;
          	//System.out.println("visited index: "+visitedIndex);
            
        	editorpane.setPage(homeURL);
        }
        if (strAction == "Close") {
            // more portable if delegated
            processWindowEvent(new WindowEvent(this,
                WindowEvent.WINDOW_CLOSING));
        }
    } catch (IOException ex) {
        ErrorMessage.showPane(null,"name: "+homeURL+"\n"+ex.getMessage(),true);
    }
}
/**
 * add buttons at the south
 */
private void addButtons() {
	btnback = new JButton("< Back");
	btnback.addActionListener(this);
	btnback.setEnabled(false);
	btnforward = new JButton("Forward >");
	btnforward.addActionListener(this);
	btnforward.setEnabled(false);
    JButton btnhome = new JButton("Home");
    btnhome.addActionListener(this);
    JButton btnclose = new JButton("Close");
    btnclose.addActionListener(this);
    //put into JPanel
    JPanel panebuttons = new JPanel();
    panebuttons.add(btnback);
    panebuttons.add(btnforward);
    panebuttons.add(btnhome);
    panebuttons.add(btnclose);
    //add panel south
    getContentPane().add(panebuttons, BorderLayout.NORTH);
}
/**
 * locate in middle of screen
 */
private void calculateLocation() {
    Dimension screendim = Toolkit.getDefaultToolkit().getScreenSize();
    setSize(new Dimension(WIDTH, HEIGHT));
    int locationx = (screendim.width - WIDTH) / 2;
    int locationy = (screendim.height - HEIGHT) / 2;
    setLocation(locationx, locationy);
}
/**
 * For testing purposes
 * @param args No argument is used.
 */
public static void main(String [] args){   
    URL index = ClassLoader.getSystemResource("index.html");
    new HelpWindow(null,"Test", index, true);

}
}