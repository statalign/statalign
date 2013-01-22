/**
 * 
 */
package statalign.ui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import statalign.base.MainManager;
import statalign.postprocess.Postprocess;

/**
 * 
 * The window for setting the input/output preferences.
 * 
 * @author miklos
 *
 */
public class OutputPreferences extends JDialog implements ActionListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	private final static Color BORDER_COLOR = Color.LIGHT_GRAY;
	
	MainFrame owner;
	
	/**
	 * Creates the window for setting the I/O preferences.
	 * 
	 * @param owner The main frame
	 */
	public OutputPreferences(MainFrame owner){
		this.owner = owner;
		Container cp = getContentPane();
		JPanel mainPanel = new JPanel(new GridLayout(1,2));
		JPanel leftPanel = new JPanel(new GridLayout(2,1));
//		leftPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(5,5,5,5),
//																BorderFactory.createEtchedBorder() 
//																/*BorderFactory.createLineBorder(Color.GRAY, 1)*/));
		mainPanel.add(leftPanel);
		JPanel rightPanel = new JPanel(new BorderLayout());
//		rightPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(5,5,5,5),
//				BorderFactory.createEtchedBorder() 
//				/*BorderFactory.createLineBorder(Color.GRAY, 1)*/));
		mainPanel.add(rightPanel);
		cp.add(mainPanel);
		//left panel
		// samplings
	//	JPanel logPanel = new JPanel(new BorderLayout());
	//	leftPanel.add(logPanel);
	//	logPanel.add(new JLabel("Select entries written into the logfile"),"North");
		List<Postprocess> postprocess = owner.manager.postProcMan.getPlugins();
		List<Postprocess> list = new ArrayList<Postprocess>();
		Comparator<Postprocess> comp = new Comparator<Postprocess>() {
			@Override
			public int compare(Postprocess o1, Postprocess o2) {
				return o1.getTabOrder()-o2.getTabOrder() < 0 ? -1 : 1;
			}
		};
		
		if(owner.manager.postProcMan.rnaMode) {
			for(Postprocess p : postprocess){
				if(p.outputable){
					list.add(p);
				}
			}
		}
		
		else {
			for(Postprocess p : postprocess) {
				if(p.outputable && !p.rnaAssociated) {
					list.add(p);
				}
			}
		}
		
		Collections.sort(list, comp);
		JPanel logListPanel = new JPanel(new GridLayout(list.size(),1));
		logListPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(5,5,5,5),
				BorderFactory.createTitledBorder(BorderFactory.createCompoundBorder(BorderFactory.createLineBorder(Color.DARK_GRAY),
						BorderFactory.createLineBorder(BORDER_COLOR)),"Select entries written into the logfile")
				/*BorderFactory.createLineBorder(Color.GRAY, 1)*/));

		leftPanel.add(logListPanel,"Center");
		JCheckBox[] logCheckBoxes = new JCheckBox[list.size()];
		int k = 0;
		for(Postprocess p : list) {
			logCheckBoxes[k] = new JCheckBox(p.getTabName(),p.sampling);
			logListPanel.add(logCheckBoxes[k]);
			logCheckBoxes[k].addActionListener(this);
			k++;
		}
		// postprocess
	//	JPanel postprocessPanel = new JPanel(new BorderLayout());
	//	leftPanel.add(postprocessPanel);
	//	postprocessPanel.add(new JLabel("Select entries that generates postprocess files"),"North");
		list.clear();
		
		if(owner.manager.postProcMan.rnaMode) {
			for(Postprocess p : postprocess){
				if(p.postprocessable){
					list.add(p);
				}
			}
		}
		
		else {
			for(Postprocess p : postprocess) {
				if(p.postprocessable && !p.rnaAssociated) {
					list.add(p);
				}
			}
		}
		
		Collections.sort(list, comp);
		JPanel postprocessListPanel = new JPanel(new GridLayout(list.size(),1));
		postprocessListPanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(5,5,5,5),
				BorderFactory.createTitledBorder(BorderFactory.createCompoundBorder(BorderFactory.createLineBorder(Color.DARK_GRAY),
						BorderFactory.createLineBorder(BORDER_COLOR)),"Select entries generating postprocess file")
				/*BorderFactory.createLineBorder(Color.GRAY, 1)*/));

		leftPanel.add(postprocessListPanel);
		JCheckBox[] postprocessCheckBoxes = new JCheckBox[list.size()];
		k = 0;
		for(Postprocess p : list) {
			postprocessCheckBoxes[k] = new JCheckBox(p.getTabName()+" ",p.postprocessWrite);
			postprocessListPanel.add(postprocessCheckBoxes[k]);
			postprocessCheckBoxes[k].addActionListener(this);
			k++;
		}
		// right panel
		//alignment types
		JRadioButton[] alignmentCheckBox = new JRadioButton[MainManager.alignmentTypes.length];
		JPanel alignmentTypePanel = new JPanel(new GridLayout(alignmentCheckBox.length,1));
		alignmentTypePanel.setBorder(BorderFactory.createCompoundBorder(BorderFactory.createEmptyBorder(5,5,5,5),
				BorderFactory.createTitledBorder(BorderFactory.createCompoundBorder(BorderFactory.createLineBorder(Color.DARK_GRAY),
						BorderFactory.createLineBorder(BORDER_COLOR)),"Select alignment output type")
				/*BorderFactory.createLineBorder(Color.GRAY, 1)*/));
		rightPanel.add(alignmentTypePanel, "Center");
		ButtonGroup alignmentGroup = new ButtonGroup();
		for(int i = 0; i < alignmentCheckBox.length; i++){
			alignmentCheckBox[i] = new JRadioButton(MainManager.alignmentTypes[i],
													owner.manager.inputData.currentAlignmentType == i);
			alignmentGroup.add(alignmentCheckBox[i]);
			alignmentCheckBox[i].addActionListener(this);
			alignmentTypePanel.add(alignmentCheckBox[i]);
		}
		
		//close button
		JPanel buttonPanel = new JPanel(new BorderLayout());
		JButton bClose = new JButton("Close");
		bClose.addActionListener(this);
		buttonPanel.add(bClose,"East");
		
		rightPanel.add(buttonPanel,"South");
		
		setTitle("Output Preferences");
		this.setBounds(owner.getX()+owner.getWidth()/10, owner.getY(), owner.getWidth()*4/5, owner.getHeight()*4/5);
		setVisible(true);
		
	}

	/* (non-Javadoc)
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
   @Override
public void actionPerformed(ActionEvent e) {
		//System.out.println(e.getSource().getClass());
		try {
			if(e.getSource().getClass() == Class.forName("javax.swing.JCheckBox")){
				//System.out.println("Found!!!");
				for(Postprocess plugin : owner.manager.postProcMan.getPlugins()){
					if(e.getActionCommand().equals(plugin.getTabName())){
						//System.out.println("This is the name: "+plugin.getTabName());
						plugin.sampling = !plugin.sampling;
					}
					if(e.getActionCommand().equals(plugin.getTabName()+" ")){
						//System.out.println("This is the name: "+plugin.getTabName());
						plugin.postprocessWrite = !plugin.postprocessWrite;
					}
				}
			}
			else if(e.getSource().getClass() == Class.forName("javax.swing.JRadioButton")){
				for(int i = 0; i < MainManager.alignmentTypes.length; i++){
					if(e.getActionCommand().equals(MainManager.alignmentTypes[i])){
						owner.manager.inputData.currentAlignmentType = i;
					}
				}
			}
			else{
				dispose();
			}
			
		} catch (ClassNotFoundException e1) {
			e1.printStackTrace();
		}
		
	}

}
