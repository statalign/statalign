package statalign.ui;

import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import statalign.base.Utils;

public class TextDlg extends JDialog {
	
	private JScrollPane scroll;

	public TextDlg(String resource, String title) throws IOException {
		BufferedReader r = new BufferedReader(new InputStreamReader(
				getClass().getClassLoader().getResource(resource).openStream()));
		String line;
		StringBuilder s = new StringBuilder();
		s.append("<html><div style='padding: 10px'><tt>");
		while((line = r.readLine()) != null) {
			s.append(line.replace(" ", "&nbsp;"));
			s.append("<br>");
		}
		s.append("</tt></div></html>");

		JEditorPane pane = new JEditorPane("text/html", s.toString());
		pane.setEditable(false);
		
		scroll = new JScrollPane(pane);
		pane.addKeyListener(new ScrollAdapter(scroll.getVerticalScrollBar(), 2, this));
		add(scroll);
		pane.setSize(new Dimension(700, 600));

		setTitle(title);
		setSize(Utils.minMax(pane.getPreferredSize().width+40, 200, 700), 600);
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
	}
	
}
