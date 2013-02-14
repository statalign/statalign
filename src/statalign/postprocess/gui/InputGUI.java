package statalign.postprocess.gui;

import java.awt.BorderLayout;
import java.awt.TextArea;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JEditorPane;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EtchedBorder;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import statalign.base.Input;
import statalign.base.MainManager;
import statalign.ui.MainFrame;

/**
 * This is the graphical interface for showing the input data
 * 
 * @author miklos, novak
 *
 */
public class InputGUI extends JPanel implements ActionListener, ListSelectionListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	Input owner;
	TextArea text;
	MainManager manager;

	JPanel seqPan;
	
	private JList sequences;
	DefaultListModel dlmSequences;
	JButton jbDelete;
	JButton jbDeleteAll;
	
	boolean showWelcome = true;
	private JScrollPane spSeq;
	
	/**
	 * This constructor makes an initial GUI for showing the input sequences and their names.
	 * @param manager The MainManager that handles the MCMC run.
	 */
	public InputGUI(final MainManager manager){
		super(new BorderLayout());
		this.manager = manager;
		
		seqPan = new JPanel(new BorderLayout());
		
		dlmSequences = new DefaultListModel();
		sequences = new JList(dlmSequences);
		sequences.setBorder(new EtchedBorder());
		sequences.setToolTipText("Input sequences - click on them to view or remove");
		sequences.addListSelectionListener(this);
		spSeq = new JScrollPane(sequences);
		//setPreferredSize(new Dimension(500,300));
//		spSzoveg.getViewport().add(sequences);
//		spSzoveg.setMaximumSize(this.getSize());
		seqPan.add(spSeq,BorderLayout.CENTER);
		
		JPanel actionPanel = new JPanel(new BorderLayout());
		jbDelete = new JButton("Remove");
		jbDelete.addActionListener(this);
		actionPanel.add(jbDelete,BorderLayout.WEST);
		
		jbDeleteAll = new JButton("Remove all");
		jbDeleteAll.addActionListener(this);
		actionPanel.add(jbDeleteAll, BorderLayout.EAST);
		
		seqPan.add(actionPanel,BorderLayout.SOUTH);

//		g.setFont(new Font("Times", Font.BOLD+Font.TRUETYPE_FONT, 16));
//		g.drawString("Welcome to StatAlign!", 30, 50);
//		g.setFont(new Font("Times", Font.BOLD, 14));
//		g.drawString("To get started, add sequences to analyse using the <b>toolbar button</b> or the File menu.", 30, 90);
//		g.drawString("If you need any more help, please refer to the documentation in the Help menu.", 30, 110);

//		JLabel welcome = new JLabel(MainFrame.WELCOME_MSG);
//		JEditorPane welcome = new JEditorPane();
//		welcome.setText(MainFrame.WELCOME_MSG);
		
		JEditorPane jep = new JEditorPane("text/html", MainFrame.WELCOME_MSG);  
		jep.setEditable(false);  
		jep.setOpaque(false);  
		jep.addHyperlinkListener(new HyperlinkListener() {  
			@Override
			public void hyperlinkUpdate(HyperlinkEvent hle) {  
				if (HyperlinkEvent.EventType.ACTIVATED.equals(hle.getEventType())) {
					String url = hle.getURL().toString();
					if(url.endsWith("add")) {
						manager.frame.addSequences();
					} else if(url.endsWith("doc"))
						manager.frame.helpUsers();
				}  
			}  
		});  
		add(jep, BorderLayout.NORTH);
		
		updateSequences();
	}
	
	/**
	 * It rereads sequences from MainManager
	 */
	public void updateSequences(){
		if(dlmSequences.size() > 0){
			dlmSequences.removeAllElements();
		}
		if(manager.inputData.seqs != null) {
			if(showWelcome && manager.inputData.seqs.size() > 0) {
				showWelcome = false;
				removeAll();
				add(seqPan);
				validate();
			}
//			System.out.println("sequences size: "+manager.seqs.sequences.size()+
//					" names size: "+manager.seqs.seqNames.size());		    
			for(int i = 0; i < manager.inputData.seqs.sequences.size(); i++){
				String s1 = manager.inputData.seqs.sequences.get(i);
				String s2 = "";
				int length = 0;
				for(int j = 0; j < s1.length(); j++){
					if(s1.charAt(j) != ' ' && s1.charAt(j) != '-'){
						s2 += s1.charAt(j);
						length++;
						if(length % 115 == 0){
							s2+="<br>";
						}
					}
				}
				dlmSequences.addElement("<html><font color=\"000099\">&gt; "+manager.inputData.seqs.seqNames.get(i)+"</font>\n<br><font face=\"Courier New\">"+s2+"</font></html>");
		    }
		}
		listListener();
	}
	
/*	
	public InputGUI(JPanel pan, Input inp, String s){
		text = new TextArea(s,pan.getHeight()/13,(int)(pan.getWidth()/6.6),TextArea.SCROLLBARS_BOTH);
		this.pan = pan;
		this.owner = inp;
	  	text.setFont(new Font("Monospaced",Font.PLAIN,10));
	  	text.setEditable(false);
	  	add(text);
	}
	
	  public void paintComponent(Graphics gr){
		  super.paintComponent(gr);
		  text.setColumns((int)(pan.getWidth()/6.6));
		  text.setRows(pan.getHeight()/13);
	  }
*/
	
	void listListener(){
		jbDeleteAll.setEnabled(false);
		
		if(sequences.getModel().getSize() != 0) {
			jbDeleteAll.setEnabled(true);
		}
		
		int index = sequences.getSelectedIndex();
		if(index == -1){
			jbDelete.setEnabled(false);
			updateUI();
		}
		else{
			jbDelete.setEnabled(true);
			updateUI();
		}
	}

	/**
	 * Handles removing sequences.
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		int index = sequences.getSelectedIndex();

		// This is apparently the best way to reverse a int[] array in Java... jeez.
		ArrayList<Integer> indices = new ArrayList<Integer>(sequences.getSelectedIndices().length);
		for (int i : sequences.getSelectedIndices()) {
			indices.add(i);
		}
		Collections.sort(indices, Collections.reverseOrder());

		// Removes the sequences from the model AND from the sequences arrays, i.e.
		// there is a one-to-one mapping between those.
		listListener();
		
		if("Remove".equals(e.getActionCommand())) {
			//listListener();
			for (int i : indices) {
				manager.inputData.seqs.seqNames.remove(i);
				manager.inputData.seqs.sequences.remove(i);
				dlmSequences.remove(i);
			}
	
			// Moves the selected index of the list.
			if(dlmSequences.getSize() != 0){
		    	if(index == dlmSequences.getSize()){
		    		index--;
		    	}
		    	sequences.setSelectedIndex(index);
		    }
		}
		
		else if("Remove all".equals(e.getActionCommand())) {
			manager.inputData.seqs.seqNames.clear();
			manager.inputData.seqs.sequences.clear();
			manager.deactivateRNA();
			dlmSequences.clear();
			jbDeleteAll.setEnabled(false);
		}
	}

	/**
	 * It invokes the list listener when a value changed.
	 */
	@Override
	public void valueChanged(ListSelectionEvent arg0) {
		listListener();
		
	}

}
