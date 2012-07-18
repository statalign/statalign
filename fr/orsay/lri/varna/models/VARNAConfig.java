/*
 VARNA is a tool for the automated drawing, visualization and annotation of the secondary structure of RNA, designed as a companion software for web servers and databases.
 Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
 electronic mail : Yann.Ponty@lri.fr
 paper mail : LRI, bat 490 Université Paris-Sud 91405 Orsay Cedex France

 This file is part of VARNA version 3.1.
 VARNA version 3.1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 VARNA version 3.1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with VARNA version 3.1.
 If not, see http://www.gnu.org/licenses.
 */
package fr.orsay.lri.varna.models;

import java.awt.Color;
import java.awt.Font;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import fr.orsay.lri.varna.models.rna.ModeleColorMap;

public class VARNAConfig implements Serializable, Cloneable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 2853916694420964233L;
	/**
	 * 
	 */
	public static final int MAJOR_VERSION = 3;
	public static final int MINOR_VERSION = 8;
	
	public static String getFullName()
	{
		return "VARNA "+MAJOR_VERSION+"."+MINOR_VERSION;
	}

	/**
	 * Enum types and internal classes
	 */

	public enum BP_STYLE implements Serializable {
		BP_STYLE_LW, BP_STYLE_SIMPLE, BP_STYLE_RNAVIZ, BP_STYLE_NONE;
		public String toString() {
			switch (this) {
			case BP_STYLE_LW:
				return "Leontis/Westhof";
			case BP_STYLE_SIMPLE:
				return "Simple (Straight line)";
			case BP_STYLE_RNAVIZ:
				return "RNAViz (Single dot)";
			case BP_STYLE_NONE:
				return "None";
			}
			return super.toString();
		}
	};

	/**
	 * Default values for config options
	 */

	public static final double MAX_ZOOM = 60;
	public static final double MIN_ZOOM = 0.5;
	public static final double DEFAULT_ZOOM = 1;
	public static final double MAX_AMOUNT = 2;
	public static final double MIN_AMOUNT = 1.01;
	public static final double DEFAULT_AMOUNT = 1.2;
	public static final double DEFAULT_BP_THICKNESS = 1.0;
	public static final double DEFAULT_DIST_NUMBERS = 3.0;

	public static final int DEFAULT_PERIOD = 10;

	public static final Color DEFAULT_TITLE_COLOR = Color.black;
	public static final Color DEFAULT_BACKBONE_COLOR = Color.DARK_GRAY.brighter();
	public static final Color DEFAULT_BOND_COLOR = Color.blue;
	public static final Color DEFAULT_SPECIAL_BASE_COLOR = Color.green.brighter();
	public static final Color DEFAULT_DASH_BASE_COLOR = Color.yellow.brighter();
	public static final double DEFAULT_BASE_OUTLINE_THICKNESS = 1.5;
	public static final Color BASE_OUTLINE_COLOR_DEFAULT = Color.DARK_GRAY.brighter();
	public static final Color BASE_INNER_COLOR_DEFAULT = new Color(242, 242,242);
	public static final Color BASE_NUMBER_COLOR_DEFAULT = Color.DARK_GRAY;
	public static final Color BASE_NAME_COLOR_DEFAULT = Color.black;
	
	public static final Color DEFAULT_HOVER_COLOR  =  new Color(230, 230,230);;

	public static final Color DEFAULT_BACKGROUND_COLOR = Color.WHITE;
	public static final Font DEFAULT_TITLE_FONT = new Font("SansSerif", Font.BOLD,18);
	public static final Font DEFAULT_BASE_FONT = new Font("SansSerif", Font.PLAIN, 18);
	public static final Font DEFAULT_NUMBERS_FONT = new Font("SansSerif",
			Font.BOLD, 18);
	public static final BP_STYLE DEFAULT_BP_STYLE = BP_STYLE.BP_STYLE_LW;

	public static final ModeleColorMap DEFAULT_COLOR_MAP = ModeleColorMap.defaultColorMap();
	public static final Color DEFAULT_COLOR_MAP_OUTLINE = Color.gray;
	public static final double DEFAULT_BP_INCREMENT = 0.65;

	public static double DEFAULT_COLOR_MAP_WIDTH = 80; 
	public static double DEFAULT_COLOR_MAP_HEIGHT = 30; 
	public static double DEFAULT_COLOR_MAP_X_OFFSET = 40; 
	public static double DEFAULT_COLOR_MAP_Y_OFFSET = 0; 
	public static int DEFAULT_COLOR_MAP_STRIPE_WIDTH = 3; 
	public static int DEFAULT_COLOR_MAP_FONT_SIZE = 20; 
	public static Color DEFAULT_COLOR_MAP_FONT_COLOR = Color.gray.darker(); 

	/**
	 * Various options.
	 */
	
	

	public double _colorMapHeight = DEFAULT_COLOR_MAP_HEIGHT; 
	public double _colorMapWidth = DEFAULT_COLOR_MAP_WIDTH; 
	public double _colorMapXOffset = DEFAULT_COLOR_MAP_X_OFFSET; 
	public double _colorMapYOffset = DEFAULT_COLOR_MAP_Y_OFFSET; 


	public BP_STYLE _mainBPStyle = DEFAULT_BP_STYLE;

	public double _zoom = DEFAULT_ZOOM;
	public double _zoomAmount = DEFAULT_AMOUNT;
	public double _bpThickness = 1.0;

	public double _baseThickness = DEFAULT_BASE_OUTLINE_THICKNESS;


	public Color _backboneColor = DEFAULT_BACKBONE_COLOR;
	public boolean _drawBackbone = true;
	
	public Color _hoverColor = DEFAULT_HOVER_COLOR;
	public Color _backgroundColor = DEFAULT_BACKGROUND_COLOR;
	public boolean _drawBackground = false;
	public Color _bondColor = DEFAULT_BOND_COLOR;
	public Color _titleColor = DEFAULT_TITLE_COLOR;
	public Color _specialBasesColor = DEFAULT_SPECIAL_BASE_COLOR;
	public Color _dashBasesColor = DEFAULT_DASH_BASE_COLOR;

	public Font _titleFont = DEFAULT_TITLE_FONT;
	public Font _numbersFont = DEFAULT_NUMBERS_FONT;
	public Font _fontBasesGeneral = DEFAULT_BASE_FONT;

	public String _title = "";

	public int _numPeriod = DEFAULT_PERIOD;

	public boolean _drawOutlineBase = true;
	public boolean _fillBase = true;
	public boolean _autoFit = true;
	public boolean _autoCenter = true;
	public boolean _modifiable = true;
	public boolean _errorsOn = false;
	public boolean _colorSpecialBases = false;
	public boolean _colorDashBases = false;
	public boolean _useBaseColorsForBPs = false;
	public boolean _drawnNonCanonicalBP = true;
	public boolean _drawnNonPlanarBP = true;
	public boolean _showWarnings = false;
	public boolean _comparisonMode = false;
	public boolean _flatExteriorLoop = false;
	
	// Relative distance between the center of a base and its number, expressed as a multiple of base radius  
	public double _distNumbers = DEFAULT_DIST_NUMBERS;
	
	public ModeleColorMap _cm = DEFAULT_COLOR_MAP;
	public String _colorMapCaption = "";
	public boolean _drawColorMap = false;

	
    public VARNAConfig clone ()
    {
        try
        {
            ByteArrayOutputStream out = new ByteArrayOutputStream ();
            ObjectOutputStream oout = new ObjectOutputStream (out);
            oout.writeObject (this);
            
            ObjectInputStream in = new ObjectInputStream (
                new ByteArrayInputStream (out.toByteArray ()));
            return (VARNAConfig)in.readObject ();
        }
        catch (Exception e)
        {
            throw new RuntimeException ("cannot clone class [" +
                this.getClass ().getName () + "] via serialization: " +
                e.toString ());
        }
    }

}
