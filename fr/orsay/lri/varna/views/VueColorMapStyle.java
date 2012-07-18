package fr.orsay.lri.varna.views;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Arrays;
import java.util.Comparator;

import javax.swing.BoxLayout;
import javax.swing.JColorChooser;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.components.GradientEditorPanel;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.rna.ModeleColorMap;
import fr.orsay.lri.varna.models.rna.ModeleColorMap.NamedColorMapTypes;

public class VueColorMapStyle extends JPanel implements ActionListener, ItemListener, PropertyChangeListener {
	private VARNAPanel _vp;
	private GradientEditorPanel _gp;
	private JComboBox _cb; 
	private JTextField _code; 
	private ModeleColorMap _backup;
	
	public VueColorMapStyle(VARNAPanel vp)
	{
		super();
		_vp = vp;
		init();
	}

	private void init()
	{
		JLabel gradientCaption = new JLabel("Click gradient to add new color...");
		_gp = new GradientEditorPanel(_vp.getColorMap().clone());
		_backup = _vp.getColorMap();
		_gp.setPreferredSize(new Dimension(300,70));
		_gp.addPropertyChangeListener(this);

		JPanel codePanel = new JPanel();
		JLabel codeCaption = new JLabel("Param. code: ");
		_code = new JTextField("");
		_code.setFont(Font.decode("Monospaced-PLAIN-12"));
		_code.setEditable(false);
		_code.addFocusListener(new FocusListener(){

			public void focusGained(FocusEvent arg0) {
						_code.setSelectionStart(0);
						_code.setSelectionEnd(_code.getText().length());
			}

			public void focusLost(FocusEvent arg0) {
			}			
		});		
		
		NamedColorMapTypes[] palettes =  ModeleColorMap.NamedColorMapTypes.values();
		Arrays.sort(palettes,new Comparator<ModeleColorMap.NamedColorMapTypes>(){
			public int compare(ModeleColorMap.NamedColorMapTypes arg0, ModeleColorMap.NamedColorMapTypes arg1) {
				return arg0.getId().compareTo(arg1.getId());
			}			
		});
		Object[] finalArray = new Object[palettes.length+1];
		int selected = -1;
		for (int i=0;i<palettes.length;i++)
		{ 
			if (palettes[i].getColorMap().equals(_vp.getColorMap()))
			{
				selected = i; 
				//System.out.println(selected);
			}
			finalArray[i] = palettes[i]; 
		}
		String custom = new String("Custom...");
		finalArray[palettes.length] = custom;
		_cb = new JComboBox(finalArray);
		if (selected!=-1)
		{
			_cb.setSelectedIndex(selected);
			_code.setText(palettes[selected].getId());
		}
		else
		{
			_cb.setSelectedItem(finalArray.length-1);
			_code.setText(_gp.getColorMap().getParamEncoding());
		}
		_cb.addItemListener(this);
		
		
		codePanel.setLayout(new BoxLayout(codePanel,BoxLayout.LINE_AXIS));
		codePanel.add(codeCaption);
		codePanel.add(_code);
		
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		add(_cb);
		add(gradientCaption);
		add(_gp);
		add(codePanel);
	}
	
	public void cancelChanges()
	{
		_vp.setColorMap(_backup);
	}
	
	public ModeleColorMap getColorMap()
	{
		return _gp.getColorMap();
	}
	
	public void actionPerformed(ActionEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	
	private void refreshCode()
	{
		int selected = -1;
		NamedColorMapTypes n = null;
		for (int i=0;i<_cb.getItemCount()-1;i++)
		{ 
			Object o = _cb.getItemAt(i);
			if (o instanceof NamedColorMapTypes)
			{
				NamedColorMapTypes ni = (NamedColorMapTypes) o;
				if (ni.getColorMap().equals(_gp.getColorMap()))
				{ 
					selected = i;	n = ni;
				}
			}
		}
		if (selected!=-1)
		{
			_code.setText(n.getId());
			_cb.setSelectedIndex(selected);
		}
		else
		{
			_code.setText(_gp.getColorMap().getParamEncoding());
		}
		_vp.setColorMap(_gp.getColorMap());
		_gp.repaint();
	}

	public void itemStateChanged(ItemEvent arg0) {
		if (arg0.getStateChange()==ItemEvent.SELECTED)
		{
		Object o = arg0.getItem();
		if (o instanceof ModeleColorMap.NamedColorMapTypes)
		{
			ModeleColorMap.NamedColorMapTypes n = ((ModeleColorMap.NamedColorMapTypes) o);
			_gp.setColorMap(n.getColorMap().clone());
			refreshCode();
		}
		}
	}

	public void propertyChange(PropertyChangeEvent arg0) {
		if (arg0.getPropertyName().equals("PaletteChanged"))
		{
			_cb.setSelectedIndex(_cb.getItemCount()-1);
			refreshCode();
		};
	}


}
