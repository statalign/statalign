/*
 VARNA is a tool for the automated drawing, visualization and annotation of the secondary structure of RNA, designed as a companion software for web servers and databases.
 Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
 electronic mail : Yann.Ponty@lri.fr
 paper mail : LRI, bat 490 Universit� Paris-Sud 91405 Orsay Cedex France

 This file is part of VARNA version 3.1.
 VARNA version 3.1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 VARNA version 3.1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with VARNA version 3.1.
 If not, see http://www.gnu.org/licenses.
 */
package fr.orsay.lri.varna.models.annotations;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;

import fr.orsay.lri.varna.models.rna.ModeleBase;
import fr.orsay.lri.varna.models.rna.VARNAPoint;



/**
 * The annotated text model
 * 
 * @author Darty@lri.fr
 * 
 */
public class TextAnnotation implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 465236085501860747L;
	/**
	 * if the annoted text is located by a static position
	 */
	public final static int POSITION = 0;
	/**
	 * if the annoted text is fixed to a base;
	 */
	public final static int BASE = 1;
	/**
	 * if the annoted text is fixed to a helix;
	 */
	public final static int HELIX = 2;
	/**
	 * if the annoted text is fixed to a loop;
	 */
	public final static int LOOP = 3;
	

	/**
	 * default text color
	 */
	public static final Color DEFAULTCOLOR = Color.black;
	/**
	 * default text font
	 */
	public static final Font DEFAULTFONT = new Font("Arial", Font.PLAIN, 12);

	private String _texte;
	private Font _font;
	private Object _ancrage;
	private int _typeAncrage;
	private Color _color;
	private double _angle;

	/**
	 * creates an annoted text on a VARNAPanel with the specified text
	 * 
	 * @param texte Textual content of the annotation
	 */
	public TextAnnotation(String texte) {
		_texte = texte;
		_color = DEFAULTCOLOR;
		_font = DEFAULTFONT;
		_angle = 0;
	}

	/**
	 * /** creates an annoted text on a VARNAPanel with the specified text and
	 * is static position
	 * 
	 * @param texte
	 * @param x
	 * @param y
	 */
	public TextAnnotation(String texte, double x, double y) {
		this(texte);
		_ancrage = new VARNAPoint(x, y);
		_typeAncrage = POSITION;
	}

	/**
	 * creates an annoted text on a VARNAPanel with the specified text fixed to
	 * a base
	 * 
	 * @param texte
	 * @param mb
	 */
	public TextAnnotation(String texte, ModeleBase mb) {
		this(texte);
		_ancrage = mb;
		_typeAncrage = BASE;
	}

	/**
	 * creates an annoted text on a VARNAPanel with the specified text fixed to
	 * a helix (if type is HELIX) or to a loop (if type is LOOP)
	 * 
	 * @param texte
	 * @param listeBase
	 * @param type
	 * @throws Exception
	 */
	public TextAnnotation(String texte, ArrayList<ModeleBase> listeBase,
			int type) throws Exception {
		this(texte);
		_ancrage = listeBase;

		if (type == HELIX)
			_typeAncrage = HELIX;
		else if (type == LOOP)
			_typeAncrage = LOOP;
		else
			throw new Exception("Bad argument");
	}

	/**
	 * creates an annoted text from another one
	 * 
	 * @param textAnnotation
	 */
	public TextAnnotation(TextAnnotation textAnnotation) {
		_ancrage = textAnnotation.getAncrage();
		_font = textAnnotation.getFont();
		_texte = textAnnotation.getTexte();
		_typeAncrage = textAnnotation.getType();
	}

	/**
	 * 
	 * @return the text
	 */
	public String getTexte() {
		return _texte;
	}

	public void setTexte(String _texte) {
		this._texte = _texte;
	}

	/**
	 * 
	 * @return the font
	 */
	public Font getFont() {
		return _font;
	}

	public void setFont(Font _font) {
		this._font = _font;
	}

	public Object getAncrage() {
		return _ancrage;
	}

	public void setAncrage(ModeleBase mb) {
		_ancrage = mb;
		_typeAncrage = BASE;
	}

	public void setAncrage(double x, double y) {
		_ancrage = new VARNAPoint(x, y);
		_typeAncrage = POSITION;
	}

	public void setAncrage(ArrayList<ModeleBase> list, int type)
			throws Exception {
		_ancrage = list;
		if (type == HELIX)
			_typeAncrage = HELIX;
		else if (type == LOOP)
			_typeAncrage = LOOP;
		else
			throw new Exception("Bad argument");
	}

	public int getType() {
		return _typeAncrage;
	}

	public Color getColor() {
		return _color;
	}

	public void setColor(Color color) {
		this._color = color;
	}

	
	public String getHelixDescription()
	{
		ArrayList<ModeleBase> listeBase =  ((ArrayList<ModeleBase>)_ancrage);
		int minA = Integer.MAX_VALUE,maxA = Integer.MIN_VALUE;
		int minB = Integer.MAX_VALUE,maxB = Integer.MIN_VALUE;
		for(ModeleBase mb : listeBase)
		{
			int i = mb.getBaseNumber();
			if (mb.getElementStructure()>i)
			{
				minA = Math.min(minA, i);
				maxA = Math.max(maxA, i);
			}
			else
			{
				minB = Math.min(minB, i);
				maxB = Math.max(maxB, i);				
			}
		}
		return "["+minA+","+maxA+"] ["+minB+","+maxB+"]";
	}
	
	public String getLoopDescription()
	{
		ArrayList<ModeleBase> listeBase =  ((ArrayList<ModeleBase>)_ancrage);
		int min = Integer.MAX_VALUE,max = Integer.MIN_VALUE;
		for(ModeleBase mb : listeBase)
		{
			int i = mb.getBaseNumber();
				min = Math.min(min, i);
				max = Math.max(max, i);
		}
		return "["+min+","+max+"]";
	}
	
	public String toString() {
		String tmp = "["+_texte+"] ";
		switch (_typeAncrage) {
		case POSITION:
			NumberFormat formatter = new DecimalFormat(".00"); 
			return tmp+" at ("+formatter.format(getCenterPosition().x)+","+formatter.format(getCenterPosition().y)+")";
		case BASE:
			return tmp+" on base "+((ModeleBase) _ancrage).getBaseNumber();
		case HELIX:
			return tmp+" on helix "+getHelixDescription();
		case LOOP:
			return tmp+" on loop "+getLoopDescription();
		default:
			return tmp;
		}		
	}

	/**
	 * 
	 * @return the text position center
	 */
	public Point2D.Double getCenterPosition() {
		switch (_typeAncrage) {
		case POSITION:
			return ((VARNAPoint) _ancrage).toPoint2D();
		case BASE:
			return ((ModeleBase) _ancrage).getCoords();
		case HELIX:
			return calculLoopHelix();
		case LOOP:
			return calculLoop();
		default:
			return new Point2D.Double(0., 0.);
		}
	}

	private Point2D.Double calculLoop() {
		ArrayList<ModeleBase> liste = extractedArrayListModeleBaseFromAncrage();
		double totalX = 0., totalY = 0.;
		for (ModeleBase base : liste) {
			totalX += base.getCoords().x;
			totalY += base.getCoords().y;
		}
		return new Point2D.Double(totalX / liste.size(), totalY / liste.size());
	}

	private Point2D.Double calculLoopHelix() {
		ArrayList<ModeleBase> liste = extractedArrayListModeleBaseFromAncrage();
		Collections.sort(liste);
		double totalX = 0., totalY = 0.;
		double num=0.0;
		for (int i=0;i<liste.size(); i++) {
			ModeleBase base =liste.get(i);
			if ((i>0 && (i<liste.size()-1)) || (liste.size()<=2))
			{
				totalX += base.getCoords().x;
				totalY += base.getCoords().y;
				num += 1;
			}
		}
		return new Point2D.Double(totalX / num, totalY / num);
	}

	
	private ArrayList<ModeleBase> extractedArrayListModeleBaseFromAncrage() {
		return (ArrayList<ModeleBase>) _ancrage;
	}

	/**
	 * clone a TextAnnotation
	 */
	public TextAnnotation clone() {
		TextAnnotation textAnnot = null;
		try {
			switch (_typeAncrage) {
			case BASE:
				textAnnot = new TextAnnotation(_texte, (ModeleBase) _ancrage);
				break;
			case POSITION:
				textAnnot = new TextAnnotation(_texte,
						((VARNAPoint) _ancrage).x,
						((VARNAPoint) _ancrage).y);
				break;
			case LOOP:
				textAnnot = new TextAnnotation(_texte,
						extractedArrayListModeleBaseFromAncrage(), LOOP);
				break;
			case HELIX:
				textAnnot = new TextAnnotation(_texte,
						extractedArrayListModeleBaseFromAncrage(), HELIX);
				break;
			default:
				break;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		textAnnot.setFont(_font);
		textAnnot.setColor(_color);
		return textAnnot;

	}

	/**
	 * copy a textAnnotation
	 * 
	 * @param textAnnotation
	 */
	public void copy(TextAnnotation textAnnotation) {
		_ancrage = textAnnotation.getAncrage();
		_font = textAnnotation.getFont();
		_texte = textAnnotation.getTexte();
		_typeAncrage = textAnnotation.getType();
		_color = textAnnotation.getColor();
		_angle = textAnnotation.getAngleInDegres();
	}


	/**
	 * 
	 * @return the angle in degrees
	 */
	public double getAngleInDegres() {
		// if (_typeAncrage == TextAnnotation.HELIX)
		// _angle = calculAngleDegres();
		return _angle;
	}

	/**
	 * 
	 * @return the angle in radians
	 */
	public double getAngleInRadians() {
		return (getAngleInDegres() * Math.PI) / 180.;
	}

	public void setAngleInDegres(double _angle) {
		this._angle = _angle;
	}

	public void setAngleInRadians(double _angle) {
		this._angle = _angle * 180 / Math.PI;
	}

}
