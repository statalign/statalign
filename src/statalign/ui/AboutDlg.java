package statalign.ui;

import java.awt.Desktop;
import java.net.URI;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

import statalign.StatAlign;
import statalign.utils.Libraries;

public class AboutDlg extends JDialog {
	
	public static final String STATAL_LICENSE = "license/STATALIGN.txt";
	public static final String VARNA_WWW = "http://";
	public static final String VARNA_LICENSE = "license/VARNA.txt";
	public static final String MPJ_WWW = "http://";
	public static final String MPJ_LICENSE = "license/MPJ.txt";
	
	private static String ABOUT_MSG;
	
	private JFrame owner;
	private JScrollPane scroll;
	private JEditorPane pane;
	
	static {
		StringBuilder msg = new StringBuilder();
		msg.append( 
			"<html>" +
		
//			"<div style='padding: 10px; background: #303030;'>" +
//			"<font face='Arial' color='white'><center>" +
////			"<img href='"+f.toURI()+"'>" +
////			"<div style='font-family: Times;font-size: 30px;color: #3526a0;'><i>StatAlign</i></div>" +
//			"<p><font size=+2>StatAlign</font><br><font size=+1>"+StatAlign.version+"</font></div>" +

			"<div style='padding: 10px'>" +
			"<font face='Arial'><center>" +
			"<p><font size=+2>StatAlign</font><br><font size=+1>"+StatAlign.version+"</font>" +
			
			"<p>(c) 2008&ndash;2013 &Aacute;d&aacute;m Nov&aacute;k et al.<br>"+
			"<a href='http://STATALIGN-www'>"+StatAlign.webPageURL+"</a><br><br>" +
			
			"<p><b>Core software by</b><br>"+
			"<p>&Aacute;d&aacute;m Nov&aacute;k<br>Istv&aacute;n Mikl&oacute;s<br>Rune Lyngs&oslash;<br>Joe Herman<br><br>" +

			"<p><b>Protein structural alignment by</b><br>"+
			"<p>Christopher Challis<br>Joe Herman<br><br>" +

			"<p><b>RNA structure prediction by</b><br>"+
			"<p>Preeti Arunapuram<br>Ing&oacute;lfur E&eth;var&eth;sson<br>Michael Golden<br><br>" +

			"<p><b>Consensus tree and<br>parallelised MC^3 by</b><br>" +
			"<p>Hrafn Eir&iacute;ksson<br>Viktor Jonsson<br>David Wood<br><br>" +

			"<p><b>Special thanks to</b><br>" +
			"<p>Jotun Hein<br>James Anderson<br>Zsuzsanna S&uuml;k&ouml;sd<br>William Majoros<br><br>" +
			
			"<p>StatAlign is distributed under the<br>" +
			"<a href='http://STATALIGN-license'>GNU General Public License Version 3</a><br><br>" +
			
			"<p>StatAlign is built using the following libraries:<br><br>" +
			"<table>");
		
		String LIBTABROW = "<tr><td>#LNAME# #VER#&ensp;<td>&ensp;<a href='http://#NAME#-www'>Web page</a>&ensp;<td>&ensp;<a href='http://#NAME#-license'>License</a>";
		for(Libraries lib : Libraries.values()) {
			if(lib.ordinal() == 0)		// skip StatAlign
				continue;
			msg.append(LIBTABROW.replaceAll("#NAME#", lib.name()).
					replaceAll("#LNAME#", lib.longName).
					replaceAll("#VER#", lib.version));
		}
		msg.append("</table></center></font></div></div></html>");
		ABOUT_MSG = msg.toString();
	}
	
	public AboutDlg(JFrame owner) {
		this.owner = owner;
		
		pane = new JEditorPane("text/html", ABOUT_MSG);  
		pane.setEditable(false);  
//		pane.setOpaque(false);  
		pane.addHyperlinkListener(new HyperlinkListener() {  
			@Override
			public void hyperlinkUpdate(HyperlinkEvent hle) {  
				if (HyperlinkEvent.EventType.ACTIVATED.equals(hle.getEventType())) {
					String url = hle.getURL().toString();
					try {
						Matcher m = Pattern.compile("http://(.*)-www").matcher(url);
						if(m.matches()) {
							URI uri = Libraries.valueOf(m.group(1)).getWebPageURI();
							Desktop.getDesktop().browse(uri);
						}
						m = Pattern.compile("http://(.*)-license").matcher(url);
						if(m.matches()) {
							Libraries lib = Libraries.valueOf(m.group(1));
							JDialog dlg = new TextDlg(lib.getLicenseFile(), lib.longName+"'s license");
							dlg.setVisible(true);
						}
						if(url.endsWith("doc")) {
						}
					} catch (Exception ex) {
						ErrorMessage.showPane(AboutDlg.this.owner, ex, true);
					}
				}  
			}  
		});
		
		scroll = new JScrollPane(pane);
		pane.addKeyListener(new ScrollAdapter(scroll.getVerticalScrollBar(), 2, this));
		add(scroll);

		setTitle("About StatAlign");
//		setModal(true);
		setSize(500, 400);
		setLocationRelativeTo(null);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
	}
	
	@Override
	public void setVisible(boolean b) {
		super.setVisible(b);
		if(b == true) {
			SwingUtilities.invokeLater(new Runnable() {
				@Override
				public void run() {
					scroll.getVerticalScrollBar().setValue(scroll.getVerticalScrollBar().getMinimum());
				}
			});
		}
		scroll.getVerticalScrollBar().setValue(scroll.getVerticalScrollBar().getMinimum());
	}
	
	public static void main(String[] args) {
		AboutDlg dlg = new AboutDlg(null);
		dlg.setVisible(true);
	}
}
