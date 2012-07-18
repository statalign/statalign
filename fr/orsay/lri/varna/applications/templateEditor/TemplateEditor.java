package fr.orsay.lri.varna.applications.templateEditor;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.ComponentOrientation;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.UIManager;
import javax.swing.filechooser.FileFilter;
import javax.swing.undo.UndoManager;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.applications.FileNameExtensionFilter;
import fr.orsay.lri.varna.exceptions.ExceptionInvalidRNATemplate;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.exceptions.ExceptionXMLGeneration;
import fr.orsay.lri.varna.models.templates.RNATemplateDrawingAlgorithmException;
import fr.orsay.lri.varna.models.templates.RNATemplateMapping;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateElement;




public class TemplateEditor extends JFrame implements KeyListener, ActionListener,DropTargetListener {

	private TemplatePanel _sk;
	private VARNAPanel _vp;
	private File currentFilePath;
	private JButton saveButton;

	public TemplateEditor()
	{
		super("VARNA Template Editor: New file");
		init();
	}
	
	private UndoManager manager;
	
	
	private void init()
	{
		try {
			_vp = new VARNAPanel(" ",".");
		} catch (ExceptionNonEqualLength e) {
			e.printStackTrace();
		}
		_vp.setNumPeriod(0);
		JPanel p = new JPanel();
		p.setLayout(new GridLayout(1,2));
		
		JToolBar systemBar = new JToolBar();
		JButton loadButton = new JButton("Open...",UIManager.getIcon("FileView.directoryIcon")); 
		loadButton.setActionCommand("open");
		loadButton.addActionListener(this);
		loadButton.addKeyListener(this);
		saveButton = new JButton("Save",UIManager.getIcon("FileView.floppyDriveIcon")); 
		saveButton.setActionCommand("save");
		saveButton.addActionListener(this);
		saveButton.addKeyListener(this);
		saveButton.setEnabled(false);
		JButton saveAsButton = new JButton("Save As...",UIManager.getIcon("FileView.floppyDriveIcon")); 
		saveAsButton.setActionCommand("save as");
		saveAsButton.addActionListener(this);
		saveAsButton.addKeyListener(this);
		JButton undoButton = new JButton("Undo"); 
		undoButton.setActionCommand("undo");
		undoButton.addActionListener(this);
		undoButton.addKeyListener(this);
		JButton redoButton = new JButton("Redo"); 
		redoButton.setActionCommand("redo");
		redoButton.addActionListener(this);
		redoButton.addKeyListener(this);
		JButton applyButtons[] = new JButton[4];
		String applyMethods[] = {"no helix length adjustment", "max increase factor", "min non-intersect factor", "helix translate"};
		for (int i=0; i<applyButtons.length; i++) {
			applyButtons[i] = new JButton("Apply (" + applyMethods[i] + ")"); 
			applyButtons[i].setActionCommand("apply" + i);
			applyButtons[i].addActionListener(this);
			applyButtons[i].addKeyListener(this);
		}
		
		JButton flipButton = new JButton("Flip"); 
		flipButton.setActionCommand("flip");
		flipButton.addActionListener(this);
		flipButton.addKeyListener(this);

		
		systemBar.add(loadButton);
		systemBar.add(saveButton);
		systemBar.add(saveAsButton);
		systemBar.addSeparator();
		systemBar.add(undoButton);
		systemBar.add(redoButton);
		systemBar.addSeparator();
		systemBar.add(flipButton);
		systemBar.addSeparator();
		for (int i=0; i<applyButtons.length; i++)
			systemBar.add(applyButtons[i]);
		systemBar.addKeyListener(this);

		JToolBar toolBar = new JToolBar();
		ButtonGroup bg = new ButtonGroup(); 
		toolBar.setOrientation(JToolBar.VERTICAL);
		//toolBar.setLayout(new BoxLayout(toolBar, BoxLayout.PAGE_AXIS));
		JToggleButton selectButton = new JToggleButton("Select"); 
		selectButton.setActionCommand("select");
		selectButton.addActionListener(this);
		selectButton.addKeyListener(this);
		JToggleButton helixButton = new JToggleButton("Helix"); 
		helixButton.setActionCommand("helix");
		helixButton.addActionListener(this);
		helixButton.addKeyListener(this);
		JToggleButton unpairedButton = new JToggleButton("Unpaired"); 
		unpairedButton.setActionCommand("unpaired");
		unpairedButton.addActionListener(this);
		unpairedButton.addKeyListener(this);

		bg.add(selectButton);
		bg.add(helixButton);
		bg.add(unpairedButton);

		toolBar.add(selectButton);
		toolBar.add(helixButton);
		toolBar.add(unpairedButton);
		systemBar.addKeyListener(this);

		
		this.setLayout(new BorderLayout());
		_sk = new TemplatePanel();
		_sk.setPreferredSize(new Dimension(800,600));
		manager = new UndoManager();
		manager.setLimit(2000);
	   _sk.addUndoableEditListener(manager);
	   _sk.addKeyListener(this);
		
		JScrollPane jp = new JScrollPane(_sk,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		p.add(jp);
		p.add(_vp);
		getContentPane().add(systemBar,BorderLayout.PAGE_START);
		getContentPane().add(toolBar,BorderLayout.WEST);
		getContentPane().add(p,BorderLayout.CENTER);
		this.addKeyListener(this);
		
	    new DropTarget(_vp, this);
	    new DropTarget(_sk, this);

		_sk.requestFocusInWindow();
	}
	
	
	private File getCurrentFilePath() {
		return currentFilePath;
	}

	private void setCurrentFilePath(File currentFilePath) {
		this.currentFilePath = currentFilePath;
		saveButton.setEnabled(true);
		setTitle("VARNA Template Editor: " + currentFilePath);
	}
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -5942729520783690050L;

	public static void main(String[] argv)
	{
		TemplateEditor frame = new TemplateEditor();
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	public void keyPressed(KeyEvent e) {
		System.out.println(e);
		switch (e.getKeyCode())
		{
		  case (KeyEvent.VK_DELETE):
		  {
			  GraphicalTemplateElement h = _sk.getSelected();
			  _sk.Unselect();
			  _sk.getTemplateUI().removeElementUI(h);
			  _sk.repaint();
		  }
		  break;
		  case (KeyEvent.VK_Z):
		  {
			  if (e.isControlDown())
			  {
				  undo();
			  }
		  }
		  break;
		  case (KeyEvent.VK_Y):
		  {
			  if (e.isControlDown())
			  {
				  redo();
			  }
		  }
		  break;
		}
	}
	
	public void undo()
	{
		if (manager.canUndo())
		{
			System.out.println("Undo: "+manager.getUndoPresentationName());
			manager.undo();
		}
	}

	public void redo()
	{
		  if (manager.canRedo())
		  {
			  System.out.println("Redo: "+manager.getRedoPresentationName());
		      manager.redo();
		  }
	}	
	
	public void keyReleased(KeyEvent e) {
		System.out.println(e);
	}

	public void keyTyped(KeyEvent e) {
		System.out.println(e);
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("undo"))
		{
			undo();
		}
		else if (e.getActionCommand().equals("redo"))
		{
			redo();
		}
		else if (e.getActionCommand().equals("flip"))
		{
			GraphicalTemplateElement gr = _sk.getSelected();
			if (gr != null)
			{
				if (gr instanceof Helix)
				{
					_sk.getTemplateUI().flipHelixUI((Helix)gr);
				}
			}
		}
		else if (e.getActionCommand().startsWith("apply"))
		{
			try {
				int method = Integer.parseInt(e.getActionCommand().substring(e.getActionCommand().length()-1));
				RNATemplateMapping map =  _vp.getRNA().drawRNATemplate(_sk.getTemplate(), method);
				for(int i: map.getSourceElemsAsSet())
				{
					RNATemplateElement t = map.getPartner(i);
					Color c= _sk.getElement(t).getDominantColor();
					_vp.getRNA().getBaseAt(i).getStyleBase().set_base_inner_color(c);
				}
				_vp.repaint();
			} catch (RNATemplateDrawingAlgorithmException e1) {
				e1.printStackTrace();
				JOptionPane.showMessageDialog(this, e1.getMessage(), "Template-based RNA drawing error", JOptionPane.ERROR_MESSAGE);
			}
		}
		else if (e.getActionCommand().equals("save"))
		{
			try {
				_sk.getTemplate().toXMLFile(getCurrentFilePath());
				System.out.println("Template saved in " + getCurrentFilePath().toString());
			} catch (ExceptionXMLGeneration e1) {
				e1.printStackTrace();
			} catch (ExceptionInvalidRNATemplate e1) {
				e1.printStackTrace();
			}
		}
		else if (e.getActionCommand().equals("save as"))
		{
			JFileChooser chooser = new JFileChooser();
			FileFilter filter = new FileNameExtensionFilter("VARNA RNA template (.xml)", "xml");
		    chooser.setFileFilter(filter);
			if (chooser.showSaveDialog(_sk) == JFileChooser.APPROVE_OPTION) {
				String path = chooser.getSelectedFile().getAbsolutePath();
				if (!path.toLowerCase().endsWith(".xml")) {
					path = path + ".xml";
				}
				try {
					_sk.getTemplate().toXMLFile(new File(path));
					System.out.println("Template saved in " + path);
					setCurrentFilePath(new File(path));
				} catch (ExceptionXMLGeneration e1) {
					e1.printStackTrace();
				} catch (ExceptionInvalidRNATemplate e1) {
					e1.printStackTrace();
				}
			}
		}
		else if (e.getActionCommand().equals("open"))
		{
			JFileChooser chooser = new JFileChooser();
		    FileFilter filter = new FileNameExtensionFilter("VARNA RNA template (.xml)", "xml");
		    chooser.setFileFilter(filter);
			if (chooser.showOpenDialog(_sk) == JFileChooser.APPROVE_OPTION) {
				File templatePath = chooser.getSelectedFile();
				_sk.loadFromXmlFile(templatePath);
				setCurrentFilePath(templatePath);
				
				// Empty the cancel history
				manager.discardAllEdits();
			}
		}
		
	}
	
	
	public void dragEnter(DropTargetDragEvent arg0) {
	}

	public void dragExit(DropTargetEvent arg0) {
	}

	public void dragOver(DropTargetDragEvent arg0) {
	}

	public void drop(DropTargetDropEvent dtde) 
	{
	  try 
	  {
	    Transferable tr = dtde.getTransferable();
	    DataFlavor[] flavors = tr.getTransferDataFlavors();
	    for (int i = 0; i < flavors.length; i++) 
	    {
	      if (flavors[i].isFlavorJavaFileListType()) 
	      {
	    	  dtde.acceptDrop(DnDConstants.ACTION_COPY_OR_MOVE);
	    	  List list = (List) tr.getTransferData(flavors[i]);
	    	  for (int j = 0; j < list.size(); j++) 
	    	  {
	    		  Object o = list.get(j);
	    		  if (dtde.getSource() instanceof DropTarget)
	    		  {
	    			  DropTarget dt = (DropTarget) dtde.getSource();
	    			  Component c = dt.getComponent();
	    			  if (c instanceof VARNAPanel)
	    			  {
	    				  VARNAPanel vp = (VARNAPanel) c;
	    				  String path = o.toString();
	    				  vp.loadFile(path,true);
	    				  vp.repaint();
	    			  }
	    			  else if (c instanceof TemplatePanel)
	    			  {
	    				  TemplatePanel sk = (TemplatePanel) c;
	    				  String path = o.toString();
	    				  sk.loadFromXmlFile(new File(path));
	    				  sk.repaint();
	    			  }
	    		  }
	    	  }
	    	  dtde.dropComplete(true);
	    	  return;
	      }
	    }
        dtde.rejectDrop();
	 } 
	 catch (Exception e) 
	 {
		 e.printStackTrace();
	     dtde.rejectDrop();
	  }
	}

	public void dropActionChanged(DropTargetDragEvent arg0) {
	}

	
}
