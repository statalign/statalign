/*
 * VARNA is a tool for the automated drawing, visualization and annotation
 * of the secondary structure of RNA, designed as a companion software for
 * web servers and databases. Copyright (C) 2008 Kevin Darty, Alain Denise
 * and Yann Ponty. electronic mail : Yann.Ponty@lri.fr paper mail : LRI, bat
 * 490 Université Paris-Sud 91405 Orsay Cedex France
 * 
 * This file is part of VARNA version 3.1. VARNA version 3.1 is free
 * software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 * 
 * VARNA version 3.1 is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with VARNA version 3.1. If not, see http://www.gnu.org/licenses.
 */
package fr.orsay.lri.varna.views;

import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;



import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.components.AnnotationTableModel;
import fr.orsay.lri.varna.controlers.ControleurTableAnnotations;

/**
 * a view for all annoted texts on the VARNAPanel
 * 
 * @author Darty@lri.fr
 * 
 */
public class VueListeAnnotations extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * if this view is for removing annoted texts
	 */
	public static final int REMOVE = 0;
	/**
	 * if this view is for editing annoted texts
	 */
	public static final int EDIT = 1;

	private VARNAPanel _vp;
	private ArrayList<Object> data;
	private JTable table;
	private int type;
	private AnnotationTableModel specialTableModel;

	/**
	 * creates the view
	 * 
	 * @param vp
	 * @param type
	 *            (REMOVE or EDIT)
	 */
	public VueListeAnnotations(VARNAPanel vp, int type) {
		super(new GridLayout(1, 0));
		this.type = type;
		_vp = vp;
		data = new ArrayList<Object>();
		data.addAll(_vp.getListeAnnotations());
		data.addAll(_vp.getRNA().getHighlightRegion());
		data.addAll(_vp.getRNA().getChemProbAnnotations());
		createView();
	}

	private void createView() {
		specialTableModel = new AnnotationTableModel(data);
		table = new JTable(specialTableModel);
		ControleurTableAnnotations ctrl = new ControleurTableAnnotations(table,
				_vp, type);
		table.addMouseListener(ctrl);
		table.addMouseMotionListener(ctrl);
		// table.setPreferredScrollableViewportSize(new Dimension(500, 100));
		// TODO: Find equivalent in JRE 1.5
		// table.setFillsViewportHeight(true);
		// Create the scroll pane and add the table to it.
		JScrollPane scrollPane = new JScrollPane(table);

		add(scrollPane);

		UIvueListeAnnotations();
	}

	/**
	 * Create the GUI and show it. For thread safety, this method should be
	 * invoked from the event-dispatching thread.
	 */
	public void UIvueListeAnnotations() {
		JComponent newContentPane = this;
		newContentPane.setOpaque(true); 
		JOptionPane.showMessageDialog(_vp, newContentPane,
				"Annotation edition", JOptionPane.PLAIN_MESSAGE);
	}

	public ArrayList<Object> getData() {
		return data;
	}

	public void setData(ArrayList<Object> data) {
		this.data = data;
	}

	public VARNAPanel get_vp() {
		return _vp;
	}

	public JTable getTable() {
		return table;
	}

	public void setTable(JTable table) {
		this.table = table;
	}

	public AnnotationTableModel getSpecialTableModel() {
		return specialTableModel;
	}

	public void setSpecialTableModel(AnnotationTableModel specialTableModel) {
		this.specialTableModel = specialTableModel;
	}
}
