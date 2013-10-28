package statalign.postprocess.gui;

import java.awt.BorderLayout;
import java.awt.TextArea;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

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
import statalign.io.DataType;
import statalign.ui.MainFrame;

/**
 * This is the graphical interface for showing the input data
 * 
 * @author miklos, novak, herman
 *
 */
public class InputGUI extends JPanel implements ActionListener, ListSelectionListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	JPanel pan;
	Input owner;
	TextArea text;
	MainManager manager;
	JPanel seqPan;
	
	private JList sequences;
	DefaultListModel dlmSequences;
	JButton jbDelete;
	JButton jbDeleteAll;
	
	public boolean sequencesAreRemovable = true;
	
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
		
//		updateSequences();
	}
	
	/**
	 * It rereads sequences from MainManager
	 */
	public void updateSequences(){
		//if(dlmSequences.size() > 0){
			dlmSequences.removeAllElements();
		//}			
		if(manager.inputData.seqs != null){
			if(showWelcome && manager.inputData.seqs.size() > 0) {
		        showWelcome = false;
		        removeAll();
		        add(seqPan);
		        validate();
			}
//			System.out.println("sequences size: "+manager.seqs.sequences.size()+
//					" names size: "+manager.seqs.seqNames.size());		    
			for(int i = 0; i < manager.inputData.seqs.size(); i++){
				String s1 = manager.inputData.seqs.getSequence(i);
				String seqTitle = "<font color=\"000099\">&gt; "+manager.inputData.seqs.getSeqName(i)+"</font>";
				for (DataType d : manager.inputData.auxData) {
					if (d.perSequenceData() && !d.getSummaryAssociatedWith(manager.inputData.seqs.getSeqName(i)).isEmpty()) {
						seqTitle += "<font color=\"C80000\"> + "+d.getSummaryAssociatedWith(manager.inputData.seqs.getSeqName(i))+"</font>";
					}
				}
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
				dlmSequences.addElement("<html>"+ seqTitle +"\n<br><font face=\"Courier New\">"+s2+"</font></html>");
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
			jbDeleteAll.setEnabled(sequencesAreRemovable);
		}
		
		int index = sequences.getSelectedIndex();
		if(index == -1){
			jbDelete.setEnabled(false);
			updateUI();
		}
		else{
			jbDelete.setEnabled(sequencesAreRemovable);
			updateUI();
		}
	}

	/**
	 * Handles removing sequences.
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		int index = sequences.getSelectedIndex();

		int[] indices = sequences.getSelectedIndices();

		// Removes the sequences from the model AND from the sequences arrays, i.e.
		// there is a one-to-one mapping between those.
		listListener();
		
		if("Remove".equals(e.getActionCommand())) {
			//listListener();
			for (int i = indices.length-1; i >= 0; i--) {
				int ind = indices[i];
				for (DataType d : manager.inputData.auxData) {
					if (d.perSequenceData()) {
						d.removeDataAssociatedWith(manager.inputData.seqs.getSeqName(ind));
					}
				}
				manager.inputData.seqs.remove(ind);
				dlmSequences.remove(ind);
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
			for (String seqName : manager.inputData.seqs.getSeqnames()) {
				for (DataType d : manager.inputData.auxData) {
					if (d.perSequenceData()) {
						d.removeDataAssociatedWith(seqName);
					}
				}
			}
			manager.inputData.seqs.clear();			
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
